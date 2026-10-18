"""P6 training: one function for every method of the common protocol (plan §4.2), driven by a method spec from
configs/p6_e1.yaml.

Every student = guided EDSR + LFC, trained N steps with L_W (Wald loss on the same patch: D(x), D(g) -> x) plus
the method's loss on its training target. Teachers (PISR / P6) are trained first on the same target and frozen;
the student is initialised from the teacher decoder when `init_from_teacher`. The final checkpoint is used
(no checkpoint selection), and evaluation uses the clean target.
"""
from __future__ import annotations

import json
import os
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader

from ..bands import GUIDE_BANDS, TARGET_BANDS
from ..degrade import band_sigmas, degrade
from . import losses as L
from . import metrics as M
from .data import SPatches, STiles
from .models import FeatureTap, Teacher, build_student


def fdl_features(net, img, guide, blocks=(2, 4)):
    """Fixed encoder for FDL on 6-band images (no VGG for 6 bands): the head and first residual blocks of a
    frozen Wald-trained student, applied to the 10 m image concatenated with the guide (its native input
    layout). Returns the features after the given numbers of blocks."""
    inner = net.model if hasattr(net, "model") else net
    h = inner.head(torch.cat([img, guide], dim=1))
    out = []
    for k, blk in enumerate(inner.body):
        h = blk(h)
        if k + 1 in blocks:
            out.append(h.float())
        if k + 1 >= max(blocks):
            break
    return out


def _save_resume(path: Path, **state):
    """Durable atomic save of the training state: write a temporary file, fsync it, rename, fsync the directory.
    Without the fsyncs a power cut after the rename can leave an empty file (seen on the GPU server, 14/10/2026)."""
    tmp = path.with_suffix(".tmp")
    with open(tmp, "wb") as f:
        torch.save(state, f)
        f.flush()
        os.fsync(f.fileno())
    tmp.replace(path)
    dfd = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(dfd)
    finally:
        os.close(dfd)


def _load_resume(path: Path, device):
    """Training state to resume from, or None. A corrupt file is moved aside and training restarts."""
    if not path.exists():
        return None
    try:
        return torch.load(path, map_location=device)
    except Exception:  # truncated by a crash: keep it for inspection, start again
        path.replace(path.with_suffix(".corrupt"))
        return None


def set_precision(cfg: dict):
    tf32 = cfg["training"].get("tf32", False)
    torch.backends.cuda.matmul.allow_tf32 = tf32
    torch.backends.cudnn.allow_tf32 = tf32


def set_seed(seed: int):
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def cosine_lr(step, total, base, warmup):
    if step < warmup:
        return base * (step + 1) / warmup
    return base * 0.5 * (1 + math.cos(math.pi * (step - warmup) / max(1, total - warmup)))


def _loader(ds, batch, workers, seed):
    return DataLoader(ds, batch_size=batch, num_workers=workers, drop_last=True,
                      generator=torch.Generator().manual_seed(seed), persistent_workers=workers > 0)


def train_teacher(target: str, cfg: dict, paths: dict, mtf: dict, out: Path, seed: int, device,
                  iterations: int | None = None) -> Path:
    """PISR teacher on the given target (privileged input): L1(recon, y_t) + lambda_T L1(E(y_t), D(y_t))."""
    out.mkdir(parents=True, exist_ok=True)
    ck = out / "teacher.pt"
    if ck.exists():
        return ck
    set_seed(seed)
    t = cfg["training"]
    teacher = Teacher(mtf).to(device)
    ds = SPatches(paths["s_train"], paths["train_tiles"], paths.get("targets"), target, t["crop"], True,
                  seed=seed + 101)
    opt = torch.optim.Adam(teacher.parameters(), lr=t["lr"])
    n_iter = iterations or t["teacher_iterations"]
    start, rs = 0, _load_resume(out / "teacher_resume.pt", device)
    if rs is not None:
        teacher.load_state_dict(rs["model"])
        opt.load_state_dict(rs["opt"])
        start = rs["step"]
    ds.offset = start * t["batch"]
    it = iter(_loader(ds, t["batch"], t["workers"], seed))
    amp = t.get("amp_bf16", False) and device.type == "cuda"
    set_precision(cfg)
    for step in range(start, n_iter):
        for pg in opt.param_groups:
            pg["lr"] = cosine_lr(step, n_iter, t["lr"], t["warmup"])
        b = {k: v.to(device, non_blocking=True) for k, v in next(it).items()}
        with torch.autocast("cuda", dtype=torch.bfloat16, enabled=amp):
            c, rec = teacher(b["target"], b["g"])
        c, rec = c.float(), rec.float()
        loss = (L.l1_masked(rec, b["target"], b["valid"])
                + t["lambda_teacher_imitation"] * (c - teacher.imitation_target(b["target"])).abs().mean())
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        if (step + 1) % t.get("resume_every", 2000) == 0:
            _save_resume(out / "teacher_resume.pt", model=teacher.state_dict(), opt=opt.state_dict(), step=step + 1)
    torch.save(teacher.state_dict(), ck)
    (out / "teacher_resume.pt").unlink(missing_ok=True)
    return ck


def train_student(method: dict, cfg: dict, paths: dict, mtf: dict, out: Path, seed: int, device,
                  iterations: int | None = None) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    if (out / "student.pt").exists() and not (out / "resume.pt").exists():
        return out / "student.pt"  # finished earlier
    set_seed(seed)
    t = cfg["training"]
    n_iter = iterations or t["iterations"]
    amp = t.get("amp_bf16", False) and device.type == "cuda"
    set_precision(cfg)
    student = build_student(mtf).to(device)
    sig, sig_g = band_sigmas(mtf, TARGET_BANDS), band_sigmas(mtf, GUIDE_BANDS)
    teacher, tap_t = None, None
    if method.get("teacher"):
        teacher = Teacher(mtf).to(device)
        teacher.load_state_dict(torch.load(paths["teacher"], map_location=device))
        teacher.eval().requires_grad_(False)
        tap_t = FeatureTap(teacher.decoder)
        if method.get("init_from_teacher", True):
            student.load_state_dict(teacher.decoder.state_dict())
    if method.get("init_ckpt"):  # e.g. B-sim+FT: second stage starts from the first-stage student
        student.load_state_dict(torch.load(method["init_ckpt"], map_location=device))
    tap_s = FeatureTap(student)
    aux = torch.nn.ModuleDict()
    w = method.get("weights", {})
    if "kd" in w:
        aux["vid"] = torch.nn.ModuleList([L.VIDEstimator(80, 80) for _ in range(2)])
    if "pcrd" in w:
        aux["ps"], aux["pt"] = L.ProjectionHead(80), L.ProjectionHead(80)
    aux = aux.to(device)
    fdl_net = None
    if "fdl" in w:
        fdl_net = build_student(mtf).to(device)
        fdl_net.load_state_dict(torch.load(paths["wald_model"], map_location=device))
        fdl_net.eval().requires_grad_(False)
    params = list(student.parameters()) + list(aux.parameters())
    opt = torch.optim.Adam(params, lr=t["lr"])
    ds = SPatches(paths["s_train"], paths["train_tiles"], paths.get("targets"), method["target"], t["crop"],
                  True, seed=seed)
    log, start = [], 0
    rs = _load_resume(out / "resume.pt", device)
    if rs is not None:  # continue an interrupted run (e.g. after a reboot)
        student.load_state_dict(rs["model"])
        aux.load_state_dict(rs["aux"])
        opt.load_state_dict(rs["opt"])
        start, log = rs["step"], rs["log"]
    ds.offset = start * t["batch"]
    it = iter(_loader(ds, t["batch"], t["workers"], seed))
    gen = torch.Generator().manual_seed(seed + start)
    t0 = time.time()
    for step in range(start, n_iter):
        for pg in opt.param_groups:
            pg["lr"] = cosine_lr(step, n_iter, t["lr"], t["warmup"])
        b = {k: v.to(device, non_blocking=True) for k, v in next(it).items()}
        x, g, yt, m = b["x"], b["g"], b["target"], b["valid"]
        terms = {}
        # Wald loss on the same patch (D(x), D(g) -> x); its features are discarded
        with torch.autocast("cuda", dtype=torch.bfloat16, enabled=amp):
            wald_pred = student(degrade(x, sig), degrade(g, sig_g))
            tap_s.pop()
            pred = student(x, g)
            fs = tap_s.pop()
        wald_pred, pred = wald_pred.float(), pred.float()
        terms["wald"] = method.get("lambda_wald", t["lambda_wald"]) * L.l1_masked(wald_pred, x)
        loss_name = method["loss"]
        if loss_name == "l1":
            terms["main"] = L.l1_masked(pred, yt, m)
        elif loss_name == "shift_l1":
            terms["main"] = L.shift_l1(pred, yt, m, radius=method.get("radius", 0.5))
        elif loss_name == "reg_hp":
            terms["main"] = L.reg_hp(pred, yt, m)
        elif loss_name == "none":
            pass
        else:
            raise ValueError(loss_name)
        if "cobi" in w:
            terms["cobi"] = w["cobi"] * L.cobi_patch(pred.float(), yt, m, w_s=method.get("cobi_ws", 0.05),
                                                     generator=gen)
        if fdl_net is not None:
            with torch.no_grad():
                f_t = fdl_features(fdl_net, yt, g)
            f_p = fdl_features(fdl_net, pred, g)
            terms["fdl"] = w["fdl"] * L.fdl(f_p, f_t, generator=gen)
        if teacher is not None:
            with torch.no_grad(), torch.autocast("cuda", dtype=torch.bfloat16, enabled=amp):
                teacher(yt, g)
                ft = tap_t.pop()
            if "kd" in w:
                hp = method.get("kd_highpass", True)
                f1 = (lambda f: L.feature_highpass(f)) if hp else (lambda f: f)  # noqa: E731
                terms["kd"] = w["kd"] * sum(est(f1(a.float()), f1(b_.float()), m)
                                            for est, a, b_ in zip(aux["vid"], fs, ft)) / len(fs)
            if "pcrd" in w:
                vs, vt = aux["ps"](fs[-1].float()), aux["pt"](ft[-1].float())
                terms["pcrd"] = w["pcrd"] * L.pcrd(vs, vt, m, mode=method["pcrd_mode"], generator=gen)
        loss = sum(terms.values())
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        if (step + 1) % t.get("resume_every", 2000) == 0 and step + 1 < n_iter:
            _save_resume(out / "resume.pt", model=student.state_dict(), aux=aux.state_dict(), opt=opt.state_dict(),
                         step=step + 1, log=log)
        if step + 1 in t.get("save_steps", []):
            torch.save(student.state_dict(), out / f"student_{step + 1}.pt")
        if (step + 1) % t["log_every"] == 0:
            log.append(dict(step=step + 1, sec=round(time.time() - t0, 1),
                            **{k: float(v.detach()) for k, v in terms.items()}))
    torch.save(student.state_dict(), out / "student.pt")
    pd.DataFrame(log).to_csv(out / "train_log.csv", index=False)
    (out / "resume.pt").unlink(missing_ok=True)
    return out / "student.pt"


@torch.no_grad()
def evaluate(ckpt: Path, paths: dict, mtf: dict, tiles: list[str], out_csv: Path, device,
             eval_h5: Path | None = None, ref_key: str = "y") -> pd.DataFrame:
    model = build_student(mtf).to(device)
    model.load_state_dict(torch.load(ckpt, map_location=device))
    model.eval()
    rows = []
    for b in DataLoader(STiles(eval_h5 or paths["s_val"], tiles, ref_key), batch_size=4, num_workers=2):
        pred = model(b["x"].to(device), b["g"].to(device)).cpu().double().numpy()
        for i, tile in enumerate(b["tile"]):
            y, v = b["y"][i].double().numpy(), b["valid"][i].numpy() > 0
            hc, hpsnr, ee = M.hp_corr(pred[i], y, v), M.hp_psnr(pred[i], y, v), M.edge_error(pred[i], y, v)
            hg = M.hp_gain(pred[i], y, v)
            for k, band in enumerate(TARGET_BANDS):
                rows.append(dict(tile=tile, site=tile.split("_")[1], band=band, hp_corr=hc[k], hp_psnr=hpsnr[k],
                                 edge_error=ee[k], hp_gain=hg[k]))
    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    return df


def save_json(path: Path, obj: dict):
    path.write_text(json.dumps(obj, indent=2, default=str))

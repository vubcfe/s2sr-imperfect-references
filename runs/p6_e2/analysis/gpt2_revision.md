# Second GPT-review round: exact site-level tests (nine test sites, seed-averaged site means)

Columns: mean difference, paired t-interval (95 %), sites + / −, exact two-sided sign-flip p (512 flips).

## H1 and H3 (reference y_c)

| test | endpoint | diff | t-interval | + / − | p | Holm p |
|---|---|---|---|---|---|---|
| H1 | hp_psnr | +4.2326 | [+3.2764, +5.1887] | 9 / 0 | 0.004 | 0.008 |
| H1 | hp_corr | +0.0611 | [+0.0482, +0.0739] | 9 / 0 | 0.004 | 0.008 |
| H1 | wald_psnr | -6.0714 | [-6.8351, -5.3077] | 0 / 9 | 0.004 | -- |
| H3 step 1 | hp_psnr | +0.3027 | [-0.3073, +0.9126] | 6 / 3 | 0.055 | 0.109 |
| H3 step 1 | hp_corr | +0.0035 | [-0.0029, +0.0099] | 7 / 2 | 0.066 | 0.109 |
| H3 step 2 | hp_psnr | +0.2755 | [-0.3510, +0.9021] | 7 / 2 | 0.355 | 0.355 |
| H3 step 2 | hp_corr | -0.0076 | [-0.0141, -0.0012] | 0 / 9 | 0.004 | 0.008 |

## Full-resolution consistency PSNR(D(y_hat), x), paired site differences (dB)

| comparison | diff | t-interval | + / − | p |
|---|---|---|---|---|
| reg_hp − reg_l1 | +1.81 | [+0.79, +2.83] | 8 / 1 | 0.012 |
| reg_hp − p6 | +2.38 | [+1.37, +3.39] | 8 / 1 | 0.008 |
| reg_hp − shift_l1 | +2.70 | [+2.02, +3.37] | 9 / 0 | 0.004 |
| reg_hp_w30 − reg_l1 | +2.69 | [+1.79, +3.58] | 9 / 0 | 0.004 |
| reg_l1 − naive_l1 | +7.18 | [+6.25, +8.10] | 9 / 0 | 0.004 |
| reg_l1 − reg_l1_nogain | +6.98 | [+6.04, +7.91] | 9 / 0 | 0.004 |

## 10 m positive control: medians of site medians (40 sites) with site-bootstrap 95 % CI

| reference | HP-gain | HP-corr | HP-PSNR | edge error |
|---|---|---|---|---|
| r | 0.616 [0.567, 0.650] | 0.778 [0.712, 0.825] | 44.804 [44.022, 45.632] | 1.397 [1.263, 1.488] |
| r_reg | 0.613 [0.558, 0.650] | 0.843 [0.820, 0.868] | 45.862 [44.720, 47.652] | 1.205 [1.105, 1.336] |
| y_c | 0.807 [0.766, 0.827] | 0.843 [0.820, 0.868] | 46.618 [45.161, 48.316] | 1.205 [1.103, 1.340] |

Compensation effect on HP-gain per site: higher at 40 of 40 sites.

Frequency separation: the 20 m high-pass passes half of the signal at 0.0094 cycles/m, where the 60 m low-pass used for a_hat transmits 0.0020.

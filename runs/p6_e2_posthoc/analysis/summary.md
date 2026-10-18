# E2 post-hoc: Wald-loss weight for Reg-HP (configs/p6_e2_posthoc.yaml; exploratory, decided after E2)

Selection on the 4 validation sites (seed 0): Reg-L1 val Wald PSNR 47.524 dB -> threshold 47.024 dB.
lambda 1: 42.36, 3: 42.98, 10: 45.21, 30: 47.31 dB -> only lambda = 30 qualifies (val HP-corr 0.8619 vs Reg-L1 0.8595).

Test sites, lambda = 30, 3 seeds (reference y): HP-PSNR 41.383 ± 0.092, HP-corr 0.8354 ± 0.0008,
HP-gain 1.127 ± 0.011, edge error 1.233 ± 0.004, Wald PSNR 47.467 ± 0.054.

Paired site bootstrap (9 sites, 2,000 resamples, 95 % CI; * = excludes 0), reference y:
| reg_hp_w30 minus | HP-PSNR | HP-corr | edge error | Wald PSNR |
|---|---|---|---|---|
| reg_l1   | -0.216 [-0.764, +0.524] | +0.0061 [+0.0010, +0.0135]* | -0.038 [-0.079, -0.012]* | -0.219 [-0.426, -0.031]* |
| shift_l1 | -0.187 [-0.687, +0.510] | +0.0085 [+0.0025, +0.0167]* | -0.035 [-0.071, -0.012]* | -0.156 [-0.381, +0.091] |
| p6       | -0.341 [-0.719, +0.036] | +0.0066 [+0.0023, +0.0116]* | -0.028 [-0.053, -0.012]* | -0.111 [-0.194, -0.014]* |
| reg_hp (lambda 1) | -0.065 [-0.212, +0.040] | -0.0010 [-0.0022, +0.0001] | +0.005 [-0.002, +0.015] | +5.470 [+4.863, +6.142]* |
| wald     | +4.167 [+3.606, +4.883]* | +0.0601 [+0.0513, +0.0708]* | -0.154 [-0.195, -0.118]* | -0.602 [-0.849, -0.431]* |

Reading: a larger Wald weight removes almost all of Reg-HP's low-frequency inconsistency (-5.5 dB -> -0.2 dB vs
Reg-L1) and keeps its structural advantage (HP-corr +0.006, edge error -0.04 vs Reg-L1); HP-PSNR not different.
Under the as-delivered reference r_reg, HP-PSNR again favours the L1 methods (-1.47 dB), the amplitude effect.
H1 (pre-specified for lambda = 1) remains "not supported"; this result is reported as post hoc.

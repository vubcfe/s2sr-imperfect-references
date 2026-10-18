# Analyses added after the GPT-5.5 review (18/10/2026)

## Full-resolution consistency: D(y_hat) vs observed 20 m bands x (test sites; mean of site means; mean and sd over seeds)

| method | PSNR (dB) | sd | contrast sd(D(y_hat))/sd(x) | corr |
|---|---|---|---|---|
| wald | 56.88 | 0.04 | 0.996 | 0.9990 |
| b_sim | 51.21 | 0.02 | 0.979 | 0.9965 |
| reg_hp_w10 | 49.22 | nan | 0.991 | 0.9946 |
| reg_hp_w30 | 48.72 | 0.13 | 0.980 | 0.9939 |
| reg_hp_w3 | 48.51 | nan | 0.988 | 0.9938 |
| reg_hp | 47.84 | 0.02 | 0.980 | 0.9928 |
| reg_l1 | 46.03 | 0.03 | 0.949 | 0.9883 |
| b_sim_ft | 46.01 | 0.15 | 0.950 | 0.9878 |
| p6 | 45.46 | 0.15 | 0.946 | 0.9863 |
| bicubic | 45.38 | nan | 0.947 | 0.9902 |
| shift_l1 | 45.14 | 0.09 | 0.949 | 0.9857 |
| reg_l1_nogain | 39.06 | 0.03 | 0.828 | 0.9597 |
| naive_l1 | 38.85 | 0.04 | 0.824 | 0.9569 |

Site-level paired differences of full-resolution PSNR (dB):

- reg_hp − reg_l1: mean +1.81, sites + / −: 8 / 1
- reg_l1 − naive_l1: mean +7.18, sites + / −: 9 / 0
- reg_l1 − reg_l1_nogain: mean +6.98, sites + / −: 9 / 0
- reg_hp_w30 − reg_l1: mean +2.69, sites + / −: 9 / 0

## Bicubic baseline (y_hat = bicubic upsampling of x), high-pass metrics

- against y: hp_psnr 40.0101, hp_corr 0.7326, edge_error 1.6925, hp_gain 0.7465
- against r_reg: hp_psnr 42.5895, hp_corr 0.7327, edge_error 1.6925, hp_gain 1.3352
- against r_raw: hp_psnr 42.2010, hp_corr 0.7026, edge_error 1.7307, hp_gain 1.3258

## Positive control on the 10 m bands (519 tiles, 40 sites): reference versions vs observed Sentinel-2 10 m (medians over tile-bands)

| reference | HP-gain | HP-corr | HP-PSNR | edge error |
|---|---|---|---|---|
| r | 0.623 | 0.827 | 44.64 | 1.272 |
| r_reg | 0.618 | 0.871 | 45.30 | 1.127 |
| y_c | 0.815 | 0.871 | 46.07 | 1.127 |

Per band (median HP-gain r -> y_c; a_hat): B02 0.61->0.78 (a_hat 0.81), B03 0.59->0.78 (a_hat 0.78), B04 0.68->0.84 (a_hat 0.80), B08 0.60->0.86 (a_hat 0.69)

Sites where y_c has higher median HP-PSNR than r_reg: 33 of 40

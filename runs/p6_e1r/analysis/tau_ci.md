# E1-R ranking reversal: Kendall tau with uncertainty

| model set | n models | endpoint | refs | tau | CI (tile bootstrap) | p (permutation) |
|---|---|---|---|---|---|---|
| all | 22 | hp_psnr | y vs y_h | 0.532 | [0.472, 0.628] | 0.0006 |
| all | 22 | hp_psnr | y vs y_inv | 1.000 | [0.991, 1.000] | 0.0000 |
| all | 22 | hp_corr | y vs y_h | 1.000 | [1.000, 1.000] | 0.0000 |
| all | 22 | hp_corr | y vs y_inv | 1.000 | [1.000, 1.000] | 0.0000 |
| one per method | 13 | hp_psnr | y vs y_h | 0.205 | [0.077, 0.333] | 0.1888 |
| one per method | 13 | hp_psnr | y vs y_inv | 1.000 | [0.974, 1.000] | 0.0000 |
| one per method | 13 | hp_corr | y vs y_h | 1.000 | [1.000, 1.000] | 0.0000 |
| one per method | 13 | hp_corr | y vs y_inv | 1.000 | [1.000, 1.000] | 0.0000 |

Ranks by HP-PSNR (tile means, seeds averaged; all models):

| model | rank vs y | rank vs y_h |
|---|---|---|
| S2_realistic/shift_l1 | 1 | 2 |
| E1R/reg_hp_h | 2 | 11 |
| S0_clean/shift_l1 | 3 | 3 |
| S2_realistic/pisr | 4 | 6 |
| S2_realistic/p6_local | 5 | 5 |
| S2_realistic/reg_l1 | 6 | 9 |
| E1R/reg_l1_clean | 7 | 8 |
| S2_realistic/p6_global | 8 | 10 |
| E1R/reg_l1_inv | 9 | 15 |
| S0_clean/pisr | 10 | 14 |
| S0_clean/p6_local | 11 | 7 |
| S0_clean/naive_l1 | 12 | 16 |
| S0_clean/reg_l1 | 13 | 13 |
| S0_clean/p6_global | 14 | 12 |
| S2_realistic/naive_l1 | 15 | 17 |
| E1R/reg_l1_h | 16 | 1 |
| S2_realistic/cobi | 17 | 19 |
| S0_clean/fdl | 18 | 4 |
| S0_clean/cobi | 19 | 20 |
| S2_realistic/fdl | 20 | 18 |
| S0_clean/wald | 21 | 22 |
| S2_realistic/wald | 22 | 21 |

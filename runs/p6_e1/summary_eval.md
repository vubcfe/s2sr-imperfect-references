# E1 results (eval sites, 181 tiles)

## Endpoints (mean over tiles; HP-corr with tile-bootstrap 95% CI)

| scenario     | method    |   hp_corr |   ci_lo |   ci_hi |   hp_psnr |   edge_error |
|:-------------|:----------|----------:|--------:|--------:|----------:|-------------:|
| S0_clean     | cobi      |    0.957  |  0.9469 |  0.9657 |   56.8266 |       0.5894 |
| S0_clean     | fdl       |    0.9537 |  0.9441 |  0.9622 |   57.0269 |       0.6352 |
| S0_clean     | naive_l1  |    0.967  |  0.9591 |  0.9737 |   58.1898 |       0.5433 |
| S0_clean     | p6_global |    0.9659 |  0.9576 |  0.9729 |   58.1716 |       0.5461 |
| S0_clean     | p6_local  |    0.9662 |  0.9581 |  0.9731 |   58.1977 |       0.5449 |
| S0_clean     | pisr      |    0.9668 |  0.9588 |  0.9736 |   58.2092 |       0.5415 |
| S0_clean     | reg_l1    |    0.9667 |  0.9588 |  0.9735 |   58.1519 |       0.5473 |
| S0_clean     | shift_l1  |    0.9685 |  0.9609 |  0.9748 |   58.4295 |       0.5434 |
| S0_clean     | wald      |    0.9141 |  0.8991 |  0.9276 |   52.3983 |       0.8331 |
| S1_shift     | cobi      |    0.9581 |  0.9481 |  0.9666 |   57.0281 |       0.5865 |
| S1_shift     | fdl       |    0.954  |  0.9447 |  0.9622 |   57.0093 |       0.6486 |
| S1_shift     | naive_l1  |    0.9646 |  0.9563 |  0.9716 |   57.9491 |       0.5684 |
| S1_shift     | p6_global |    0.9673 |  0.9594 |  0.974  |   58.3239 |       0.5442 |
| S1_shift     | p6_local  |    0.9674 |  0.9596 |  0.9741 |   58.317  |       0.5431 |
| S1_shift     | pisr      |    0.9682 |  0.9605 |  0.9747 |   58.4074 |       0.5417 |
| S1_shift     | reg_l1    |    0.9687 |  0.9612 |  0.975  |   58.4187 |       0.5436 |
| S1_shift     | shift_l1  |    0.9694 |  0.9622 |  0.9755 |   58.6292 |       0.5417 |
| S1_shift     | wald      |    0.914  |  0.899  |  0.9274 |   52.403  |       0.834  |
| S2_realistic | cobi      |    0.958  |  0.9481 |  0.9665 |   57.182  |       0.5879 |
| S2_realistic | fdl       |    0.9524 |  0.9426 |  0.9611 |   56.7373 |       0.644  |
| S2_realistic | naive_l1  |    0.9644 |  0.9562 |  0.9715 |   57.8768 |       0.5699 |
| S2_realistic | p6_global |    0.9669 |  0.9589 |  0.9737 |   58.2841 |       0.5472 |
| S2_realistic | p6_local  |    0.9672 |  0.9593 |  0.974  |   58.3388 |       0.5465 |
| S2_realistic | pisr      |    0.9679 |  0.9602 |  0.9745 |   58.3724 |       0.5419 |
| S2_realistic | reg_l1    |    0.9679 |  0.9602 |  0.9744 |   58.2988 |       0.5467 |
| S2_realistic | shift_l1  |    0.9692 |  0.9619 |  0.9753 |   58.5452 |       0.5418 |
| S2_realistic | wald      |    0.9138 |  0.8989 |  0.9273 |   52.3905 |       0.8339 |
| S3_stress    | cobi      |    0.9592 |  0.9496 |  0.9675 |   57.479  |       0.5832 |
| S3_stress    | fdl       |    0.9439 |  0.9329 |  0.9538 |   55.8465 |       0.6972 |
| S3_stress    | naive_l1  |    0.9587 |  0.9496 |  0.9666 |   57.4937 |       0.6122 |
| S3_stress    | p6_global |    0.9627 |  0.9537 |  0.9704 |   57.7325 |       0.5721 |
| S3_stress    | p6_local  |    0.9631 |  0.9542 |  0.9707 |   57.8309 |       0.5718 |
| S3_stress    | pisr      |    0.9633 |  0.9547 |  0.9708 |   57.7716 |       0.5702 |
| S3_stress    | reg_l1    |    0.9647 |  0.9564 |  0.9718 |   57.9306 |       0.5691 |
| S3_stress    | shift_l1  |    0.9662 |  0.9583 |  0.973  |   58.1088 |       0.5604 |
| S3_stress    | wald      |    0.9141 |  0.8992 |  0.9276 |   52.4088 |       0.8343 |

## Drop of HP-corr vs the method's own S0_clean (H2)

| scenario     | method    |    drop |   ci_lo |   ci_hi | robust   | degrades   |
|:-------------|:----------|--------:|--------:|--------:|:---------|:-----------|
| S1_shift     | cobi      | -0.0011 | -0.0014 | -0.0009 | True     | False      |
| S1_shift     | fdl       | -0.0003 | -0.0006 |  0.0001 | True     | False      |
| S1_shift     | naive_l1  |  0.0024 |  0.002  |  0.0029 | True     | False      |
| S1_shift     | p6_global | -0.0014 | -0.0018 | -0.001  | True     | False      |
| S1_shift     | p6_local  | -0.0013 | -0.0016 | -0.001  | True     | False      |
| S1_shift     | pisr      | -0.0014 | -0.0017 | -0.001  | True     | False      |
| S1_shift     | reg_l1    | -0.002  | -0.0025 | -0.0015 | True     | False      |
| S1_shift     | shift_l1  | -0.001  | -0.0013 | -0.0007 | True     | False      |
| S1_shift     | wald      |  0.0001 |  0.0001 |  0.0002 | True     | False      |
| S2_realistic | cobi      | -0.0011 | -0.0013 | -0.0008 | True     | False      |
| S2_realistic | fdl       |  0.0013 |  0.001  |  0.0016 | True     | False      |
| S2_realistic | naive_l1  |  0.0025 |  0.0022 |  0.0029 | True     | False      |
| S2_realistic | p6_global | -0.001  | -0.0013 | -0.0007 | True     | False      |
| S2_realistic | p6_local  | -0.001  | -0.0013 | -0.0007 | True     | False      |
| S2_realistic | pisr      | -0.0011 | -0.0014 | -0.0008 | True     | False      |
| S2_realistic | reg_l1    | -0.0012 | -0.0015 | -0.0009 | True     | False      |
| S2_realistic | shift_l1  | -0.0008 | -0.0011 | -0.0005 | True     | False      |
| S2_realistic | wald      |  0.0003 |  0.0002 |  0.0003 | True     | False      |
| S3_stress    | cobi      | -0.0022 | -0.0027 | -0.0017 | True     | False      |
| S3_stress    | fdl       |  0.0099 |  0.0083 |  0.0115 | True     | False      |
| S3_stress    | naive_l1  |  0.0082 |  0.007  |  0.0096 | True     | False      |
| S3_stress    | p6_global |  0.0032 |  0.0026 |  0.004  | True     | False      |
| S3_stress    | p6_local  |  0.0031 |  0.0024 |  0.0038 | True     | False      |
| S3_stress    | pisr      |  0.0035 |  0.0027 |  0.0043 | True     | False      |
| S3_stress    | reg_l1    |  0.002  |  0.0016 |  0.0025 | True     | False      |
| S3_stress    | shift_l1  |  0.0022 |  0.0018 |  0.0027 | True     | False      |
| S3_stress    | wald      | -0.0001 | -0.0001 | -0      | True     | False      |

## E2 loss selection on S2_realistic

Edge-error constraint: <= Reg-L1 + 0.05 px = 0.597; eligible: ['naive_l1', 'reg_l1', 'shift_l1', 'cobi', 'pisr', 'p6_local', 'p6_global']
**Selected: shift_l1**

| vs        |   diff |   p |   p_holm |
|:----------|-------:|----:|---------:|
| naive_l1  | 0.0048 |   0 |        0 |
| reg_l1    | 0.0013 |   0 |        0 |
| cobi      | 0.0112 |   0 |        0 |
| fdl       | 0.0168 |   0 |        0 |
| pisr      | 0.0013 |   0 |        0 |
| p6_local  | 0.002  |   0 |        0 |
| p6_global | 0.0024 |   0 |        0 |

pCRD-L vs Shift-L1: diff -0.0020, p 0.0000 -> pCRD-L does not enter E2
## Seed replication (mean HP-corr per seed; mean, sd over seeds)

| scenario     | method    |    mean |     std |   count |
|:-------------|:----------|--------:|--------:|--------:|
| S0_clean     | naive_l1  | 0.96736 | 0.0004  |       3 |
| S0_clean     | p6_global | 0.96564 | 0.00025 |       3 |
| S0_clean     | reg_l1    | 0.96712 | 0.00032 |       3 |
| S0_clean     | shift_l1  | 0.96867 | 0.00019 |       3 |
| S2_realistic | naive_l1  | 0.96479 | 0.00043 |       3 |
| S2_realistic | p6_global | 0.96664 | 0.00027 |       3 |
| S2_realistic | reg_l1    | 0.96791 | 0.00013 |       3 |
| S2_realistic | shift_l1  | 0.96959 | 0.00038 |       3 |

Paired differences across seeds (same seed, same scenario):

| scenario     | diff                |   n_seeds |     mean |      sd |      min |      max |
|:-------------|:--------------------|----------:|---------:|--------:|---------:|---------:|
| S0_clean     | shift_l1 - reg_l1   |         3 |  0.00155 | 0.00013 |  0.00146 |  0.0017  |
| S0_clean     | reg_l1 - naive_l1   |         3 | -0.00025 | 0.00024 | -0.00049 | -2e-05   |
| S0_clean     | p6_global - reg_l1  |         3 | -0.00147 | 0.00056 | -0.00189 | -0.00084 |
| S0_clean     | shift_l1 - naive_l1 |         3 |  0.00131 | 0.00029 |  0.00097 |  0.00148 |
| S2_realistic | shift_l1 - reg_l1   |         3 |  0.00168 | 0.00049 |  0.00128 |  0.00223 |
| S2_realistic | reg_l1 - naive_l1   |         3 |  0.00312 | 0.00054 |  0.0025  |  0.00351 |
| S2_realistic | p6_global - reg_l1  |         3 | -0.00126 | 0.00036 | -0.00168 | -0.00105 |
| S2_realistic | shift_l1 - naive_l1 |         3 |  0.0048  | 9e-05   |  0.00472 |  0.0049  |

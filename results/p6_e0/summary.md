# E0 summary: 2100 reliable tiles of 2269 (status ok: 2110, of which g_peak >= 0.1: 2100), 49 sites, 116 flights among reliable tiles (119 in all tiles)

Filter scales: low-pass sigma 60 m, high-pass removes sigma 20 m; blocks 640 m.

## Tile status

| status                        |   count |
|:------------------------------|--------:|
| ok                            |    2110 |
| skipped: 20 m bands resampled |     131 |
| skipped: <50% valid           |      28 |

## Geometry (10 m pixels)

g/t: guide (10 m) / target (20 m) bands; gt_diff: difference of their global shifts;
rms_translation / rms_affine: block shifts left after one translation / an affine field.

|                   |   q10 |   q25 |   q50 |   q75 |   q90 |
|:------------------|------:|------:|------:|------:|------:|
| g_mag             | 0.114 | 0.175 | 0.279 | 0.532 | 0.767 |
| t_mag             | 0.128 | 0.2   | 0.322 | 0.569 | 0.84  |
| gt_diff           | 0.03  | 0.054 | 0.095 | 0.15  | 0.261 |
| rms_translation_g | 0.031 | 0.045 | 0.076 | 0.132 | 0.221 |
| rms_affine_g      | 0.024 | 0.035 | 0.061 | 0.103 | 0.168 |
| rms_translation_t | 0.079 | 0.108 | 0.168 | 0.26  | 0.43  |
| rms_affine_t      | 0.066 | 0.094 | 0.14  | 0.214 | 0.344 |

- tiles with |global shift| > 0.25 px: 57.7%
- tiles with |global shift| > 0.5 px: 27.2%
- tiles with |global shift| > 1.0 px: 4.5%

### H0 (one translation per tile leaves <= 0.25 px)

Pre-specified measure (rms of block shifts around their median):

- guide bands (10 m): available for 2076 of 2100 tiles; blocks per tile (median, same tiles): 16; median 0.076, q90 0.221; share <= 0.25: **91.7%** (site-bootstrap CI 88.5%-94.2%); missing counted as failing: 90.6% (CI 87.4%-93.1%)
- target bands (20 m): available for 2060 of 2100 tiles; blocks per tile (median, same tiles): 15; median 0.168, q90 0.430; share <= 0.25: **73.0%** (site-bootstrap CI 64.8%-79.8%); missing counted as failing: 71.6% (CI 63.4%-78.3%)

After the shifts the correction uses (guide bands: guide tile shift; target bands: joint target-band shift); rms of block shifts around 0:

- guide bands: available for 2079 of 2100 tiles; blocks per tile (median, same tiles): 15; median 0.083, q90 0.235; share <= 0.25: **91.2%** (site-bootstrap CI 87.9%-93.9%); missing counted as failing: 90.3% (CI 86.8%-93.0%)
- target bands (all six, raw; includes estimation noise): available for 2074 of 2100 tiles; blocks per tile (median, same tiles): 14; median 0.178, q90 0.452; share <= 0.25: **71.0%** (site-bootstrap CI 62.8%-78.2%); missing counted as failing: 70.1% (CI 61.8%-77.3%)
- target group red edge B5-B8A (raw): available for 2067 of 2100 tiles; median 0.232, q90 0.771; share <= 0.25: **54.0%** (site-bootstrap CI 43.2%-64.1%); missing counted as failing: 53.2% (CI 42.2%-62.9%)
- target group SWIR B11-B12 (raw): available for 2027 of 2100 tiles; median 0.247, q90 0.595; share <= 0.25: **50.9%** (site-bootstrap CI 39.8%-61.2%); missing counted as failing: 49.1% (CI 38.1%-58.8%)

H0' (post hoc): field common to the two target groups, unbiased squared estimate common2 = mean(a . b) thresholded at 0.25^2 (the clipped sqrt `common` is biased up by ~0.04 px at 14 blocks):

- H0' common target field (common2 <= 0.0625): available for 1739 of 2100 tiles; blocks per tile (median, same tiles): 14; median 0.004, q90 0.082; share <= 0.0625: **87.2%** (site-bootstrap CI 82.4%-90.6%); missing counted as failing: 72.2% (CI 66.4%-76.9%)
- field common to target and guide (common2 <= 0.0625): available for 1739 of 2100 tiles; blocks per tile (median, same tiles): 14; median 0.005, q90 0.047; share <= 0.0625: **92.5%** (site-bootstrap CI 90.0%-94.5%); missing counted as failing: 76.6% (CI 72.3%-80.2%)
- estimation noise of one target group (rms): available for 1739 of 2100 tiles; blocks per tile (median, same tiles): 14; median 0.226, q90 0.745; share <= 0.25: **57.5%** (site-bootstrap CI 45.7%-67.8%); missing counted as failing: 47.6% (CI 37.9%-55.9%)

Sensitivity: target bands aligned with the guide shift instead:

- target bands (raw): available for 2075 of 2100 tiles; median 0.202, q90 0.583; share <= 0.25: **63.3%** (site-bootstrap CI 53.5%-72.2%); missing counted as failing: 62.5% (CI 52.8%-71.3%)
- H0' measure: available for 1739 of 2100 tiles; median 0.012, q90 0.125; share <= 0.0625: **81.9%** (site-bootstrap CI 75.1%-87.0%); missing counted as failing: 67.8% (CI 60.6%-73.5%)

Sensitivity, excluding the top cloud-score quartile (230 tiles):

- guide bands: available for 1855 of 1870 tiles; median 0.078, q90 0.208; share <= 0.25: **92.8%** (site-bootstrap CI 90.3%-95.0%); missing counted as failing: 92.1% (CI 89.5%-94.2%)
- H0' measure: available for 1570 of 1870 tiles; median 0.003, q90 0.066; share <= 0.0625: **89.2%** (site-bootstrap CI 85.3%-91.9%); missing counted as failing: 74.9% (CI 70.2%-78.4%)

Geometry by |dt| (days): if residual fields were surface change they would grow with |dt| (R3-1):

| dt_bin   |    n |   rms_tileshift_g |   rms_tileshift_t |   common_t |   H0prime_share |   common_tg |
|:---------|-----:|------------------:|------------------:|-----------:|----------------:|------------:|
| 0-3      |  220 |             0.049 |             0.105 |      0.032 |           0.973 |       0.039 |
| 4-7      |  164 |             0.069 |             0.107 |      0     |           0.979 |       0.044 |
| 8-15     |  220 |             0.08  |             0.185 |      0.096 |           0.87  |       0.072 |
| 16-40    | 1496 |             0.092 |             0.198 |      0.073 |           0.844 |       0.077 |

Edge vs interior blocks after the tile shift (medians, 10 m px; reflect-padded tile edges, R3-3.NEW-3): rms_edge_g 0.085, rms_inner_g 0.065, rms_edge_t 0.179, rms_inner_t 0.146


Phase-correlation peaks:

|        |   q10 |   q25 |   q50 |   q75 |   q90 |
|:-------|------:|------:|------:|------:|------:|
| g_peak | 0.47  | 0.618 | 0.734 | 0.814 | 0.864 |
| t_peak | 0.455 | 0.605 | 0.72  | 0.804 | 0.866 |

### Where shifts vary (share of sum of squares; R1-1.M2)

|      |   site |   flight_in_site |   tile_in_flight |
|:-----|-------:|-----------------:|-----------------:|
| g_dy |  0.225 |            0.362 |            0.413 |
| g_dx |  0.289 |            0.279 |            0.432 |

## Radiometry per band, cloud-clean subset (1870 tiles, cloud score <= q75; medians; reflectance units)

Statistics use only pixels whose whole 4-sigma low-pass footprint is valid. hp_gain is an *apparent* gain: S2 noise and aliasing in x lower it; not yet corrected (pre-E2 gate).

rmse_al: after global shift; rmse_lin_*: after shift + gain/offset fitted on the low-pass (ols, symmetric, or cross-fitted on the other tile half); hp_gain: std / OLS ratio of high-pass detail r vs x.

| band   |   rmse_raw |   rmse_al |   rmse_lin_ols |   rmse_lin_sym |   rmse_lin_cv |   a_ols |   a_sd |   r2_lowpass |   hp_gain_sd |   hp_gain_ols |   hp_corr_raw |   hp_corr |
|:-------|-----------:|----------:|---------------:|---------------:|--------------:|--------:|-------:|-------------:|-------------:|--------------:|--------------:|----------:|
| B02    |     0.0076 |    0.0073 |         0.0048 |         0.0051 |        0.0053 |  0.6781 | 0.764  |       0.8358 |       0.583  |        0.4719 |        0.7761 |    0.829  |
| B03    |     0.0096 |    0.0092 |         0.0057 |         0.0061 |        0.0063 |  0.6406 | 0.7357 |       0.8163 |       0.5634 |        0.4573 |        0.7618 |    0.8241 |
| B04    |     0.0128 |    0.0121 |         0.0077 |         0.0082 |        0.0086 |  0.6914 | 0.7582 |       0.859  |       0.6344 |        0.5389 |        0.8049 |    0.8597 |
| B08    |     0.0297 |    0.0285 |         0.0148 |         0.0158 |        0.0163 |  0.5265 | 0.6276 |       0.7575 |       0.5209 |        0.3938 |        0.707  |    0.7809 |
| B05    |     0.0136 |    0.0134 |         0.0063 |         0.0066 |        0.007  |  0.5329 | 0.6055 |       0.8117 |       0.4871 |        0.3968 |        0.792  |    0.8235 |
| B06    |     0.0233 |    0.0231 |         0.0092 |         0.0098 |        0.0103 |  0.429  | 0.5217 |       0.7448 |       0.4253 |        0.3114 |        0.7209 |    0.7577 |
| B07    |     0.027  |    0.0267 |         0.0106 |         0.0113 |        0.0122 |  0.4397 | 0.5365 |       0.749  |       0.4355 |        0.3244 |        0.7232 |    0.7599 |
| B8A    |     0.0277 |    0.0272 |         0.011  |         0.0116 |        0.0126 |  0.4391 | 0.5345 |       0.7562 |       0.4459 |        0.3343 |        0.7239 |    0.7656 |
| B11    |     0.0213 |    0.021  |         0.0091 |         0.0096 |        0.0106 |  0.5129 | 0.5763 |       0.8314 |       0.5591 |        0.4644 |        0.8073 |    0.8449 |
| B12    |     0.0177 |    0.0175 |         0.008  |         0.0084 |        0.0094 |  0.5573 | 0.6187 |       0.8474 |       0.5714 |        0.4831 |        0.8323 |    0.8609 |

RMSE reduction computed from the medians above (R1-3.7):

| band   |   align_pct (ratio of medians) |   lin_cv_pct (ratio of medians) |
|:-------|-------------------------------:|--------------------------------:|
| B02    |                            3.9 |                            27.6 |
| B03    |                            3.9 |                            31.7 |
| B04    |                            5.1 |                            28.9 |
| B08    |                            4   |                            42.8 |
| B05    |                            1.2 |                            47.6 |
| B06    |                            1.2 |                            55.4 |
| B07    |                            1.1 |                            54.2 |
| B8A    |                            1.5 |                            53.9 |
| B11    |                            1.4 |                            49.9 |
| B12    |                            1.2 |                            46.1 |

All reliable tiles incl. cloudy ones (2100; sensitivity):

| band   |   rmse_lin_cv |   a_sd |   hp_gain_sd |   hp_corr |
|:-------|--------------:|-------:|-------------:|----------:|
| B02    |        0.0053 | 0.7538 |       0.5757 |    0.8199 |
| B03    |        0.0062 | 0.7238 |       0.5552 |    0.8154 |
| B04    |        0.0085 | 0.7501 |       0.6272 |    0.8495 |
| B08    |        0.0166 | 0.612  |       0.5067 |    0.7683 |
| B05    |        0.007  | 0.5994 |       0.4797 |    0.8124 |
| B06    |        0.0104 | 0.5089 |       0.4126 |    0.7444 |
| B07    |        0.0123 | 0.5181 |       0.4221 |    0.7491 |
| B8A    |        0.0127 | 0.5166 |       0.4336 |    0.7504 |
| B11    |        0.0105 | 0.5646 |       0.5488 |    0.8362 |
| B12    |        0.0093 | 0.6088 |       0.5613 |    0.8539 |

## Harmonisation check

A per-tile regression of S2 on the reference would match tile means exactly and shrink the reference's standard deviation to about R times that of S2 at the fitted scale.

| band   |   median R |   median sd ratio |   median hp_gain_sd |   corr(sd ratio, R) over tiles |   median |mean_r - mean_x| |
|:-------|-----------:|------------------:|--------------------:|-------------------------------:|---------------------------:|
| B02    |     0.881  |            0.7164 |              0.583  |                         0.7505 |                     0.0004 |
| B03    |     0.869  |            0.6934 |              0.5634 |                         0.7272 |                     0.0006 |
| B04    |     0.9029 |            0.7395 |              0.6344 |                         0.7368 |                     0.0008 |
| B08    |     0.83   |            0.6037 |              0.5209 |                         0.729  |                     0.0018 |
| B05    |     0.8781 |            0.5851 |              0.4871 |                         0.6375 |                     0.0011 |
| B06    |     0.8246 |            0.5016 |              0.4253 |                         0.6479 |                     0.0018 |
| B07    |     0.8347 |            0.5165 |              0.4355 |                         0.6745 |                     0.0021 |
| B8A    |     0.8366 |            0.5172 |              0.4459 |                         0.6557 |                     0.0021 |
| B11    |     0.8879 |            0.5767 |              0.5591 |                         0.5889 |                     0.0019 |
| B12    |     0.9037 |            0.6148 |              0.5714 |                         0.6195 |                     0.0015 |

## Target bands: paired gain of each band's own shift over the joint target-band shift

d_rmse = rmse(tile shift) - rmse(own shift), d_hpcorr = hp_corr(own) - hp_corr(tile); positive = own shift better. Rows: tiles where the band's own peak >= min_peak.

| band   |    n |   own_shift_q50 |   own_shift_q90 |   d_rmse_median |   d_rmse_q90 |   share_rmse_gain_gt_0_001 |   d_hpcorr_median |   d_hpcorr_q90 |   share_hpcorr_gain_gt_0_01 |
|:-------|-----:|----------------:|----------------:|----------------:|-------------:|---------------------------:|------------------:|---------------:|----------------------------:|
| B05    | 1863 |          0.2953 |          0.7955 |               0 |       0.0001 |                     0.0005 |            0.0005 |         0.0054 |                      0.0424 |
| B06    | 1861 |          0.3423 |          0.9062 |               0 |       0.0001 |                     0      |            0.0007 |         0.0064 |                      0.0553 |
| B07    | 1858 |          0.3661 |          0.933  |               0 |       0.0001 |                     0      |            0.0009 |         0.009  |                      0.0861 |
| B8A    | 1859 |          0.3758 |          0.9273 |               0 |       0.0001 |                     0      |            0.0012 |         0.0095 |                      0.0898 |
| B11    | 1864 |          0.3677 |          0.8914 |               0 |       0.0001 |                     0      |            0.0027 |         0.0096 |                      0.0928 |
| B12    | 1864 |          0.3256 |          0.8796 |               0 |       0.0001 |                     0.0016 |            0.0017 |         0.0095 |                      0.0923 |

By quartile of the guide-target shift difference gt_diff (10 m px):

| gt_diff        |    n |   d_rmse_median |   d_hpcorr_median |   share_hpcorr_gain_gt_0_01 |
|:---------------|-----:|----------------:|------------------:|----------------------------:|
| (-0.001, 0.05] | 2812 |               0 |            0.0006 |                      0.0569 |
| (0.05, 0.0906] | 2782 |               0 |            0.0009 |                      0.0575 |
| (0.0906, 0.14] | 2790 |               0 |            0.0014 |                      0.0724 |
| (0.14, 3.662]  | 2785 |               0 |            0.002  |                      0.1192 |

## Cross-fitted residual RMSE by |dt| (days)

tiles per bin: 0-3: 220, 4-7: 164, 8-15: 220, 16-40: 1496

| band   |    0-3 |    4-7 |   8-15 |   16-40 |
|:-------|-------:|-------:|-------:|--------:|
| B02    | 0.0046 | 0.0052 | 0.0045 |  0.0056 |
| B03    | 0.0053 | 0.0062 | 0.0054 |  0.0066 |
| B04    | 0.0065 | 0.0083 | 0.0069 |  0.0094 |
| B08    | 0.0122 | 0.0163 | 0.017  |  0.0166 |
| B05    | 0.0052 | 0.0072 | 0.0062 |  0.0075 |
| B06    | 0.0081 | 0.0108 | 0.0104 |  0.0105 |
| B07    | 0.0093 | 0.0125 | 0.0127 |  0.0125 |
| B8A    | 0.0095 | 0.0132 | 0.0133 |  0.0128 |
| B11    | 0.0078 | 0.0102 | 0.0098 |  0.0113 |
| B12    | 0.0071 | 0.0099 | 0.0077 |  0.01   |

## Same, by land-cover superclass (medians; small classes are indicative only)

tiles per class: Rural: 1088, Forest: 757, Water: 101, Developed: 30

| band   |   Developed |   Forest |   Rural |   Water |
|:-------|------------:|---------:|--------:|--------:|
| B02    |      0.0131 |   0.0042 |  0.0058 |  0.0049 |
| B03    |      0.0142 |   0.0053 |  0.0069 |  0.0062 |
| B04    |      0.016  |   0.0063 |  0.0102 |  0.0072 |
| B08    |      0.0243 |   0.0174 |  0.0141 |  0.0304 |
| B05    |      0.0121 |   0.0057 |  0.008  |  0.007  |
| B06    |      0.0124 |   0.0099 |  0.01   |  0.0194 |
| B07    |      0.0166 |   0.0121 |  0.0112 |  0.0237 |
| B8A    |      0.0161 |   0.0125 |  0.0114 |  0.0244 |
| B11    |      0.0155 |   0.0096 |  0.0115 |  0.0132 |
| B12    |      0.015  |   0.0077 |  0.0113 |  0.0093 |

## High-pass gain (sd ratio) of target bands by land cover

| band   |   Developed |   Forest |   Rural |   Water |
|:-------|------------:|---------:|--------:|--------:|
| B05    |       0.568 |    0.441 |   0.507 |   0.506 |
| B06    |       0.426 |    0.362 |   0.441 |   0.554 |
| B07    |       0.427 |    0.384 |   0.448 |   0.585 |
| B8A    |       0.438 |    0.397 |   0.457 |   0.611 |
| B11    |       0.594 |    0.524 |   0.565 |   0.685 |
| B12    |       0.592 |    0.551 |   0.571 |   0.651 |

## By tile cloud score (cloud_score_plus_cdf quartiles; SEN2NEON has no per-pixel cloud mask). Radiometry columns: medians over the six target bands and tiles.

| cloud_score    |   n |   g_peak |   rms_tileshift_g |   rmse_lin_cv |     a_sd |   hp_gain_sd |   hp_corr |
|:---------------|----:|---------:|------------------:|--------------:|---------:|-------------:|----------:|
| (0.862, 0.906] | 540 |   0.7731 |            0.0496 |        0.0079 |   0.5497 |       0.4563 |    0.8148 |
| (0.906, 0.914] | 611 |   0.7704 |            0.0777 |        0.011  |   0.5876 |       0.5065 |    0.8169 |
| (0.914, 0.922] | 719 |   0.7054 |            0.1017 |      nan      | nan      |     nan      |  nan      |
| (0.922, 0.937] | 230 |   0.521  |            0.1348 |      nan      | nan      |     nan      |  nan      |

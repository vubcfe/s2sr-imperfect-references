# E0 summary: 419 reliable tiles of 454 (status ok: 422, of which g_peak >= 0.1: 419), 44 sites, 102 flights among reliable tiles (108 in all tiles)

Filter scales: low-pass sigma 60 m, high-pass removes sigma 20 m; blocks 640 m.

## Tile status

| status                        |   count |
|:------------------------------|--------:|
| ok                            |     422 |
| skipped: 20 m bands resampled |      26 |
| skipped: <50% valid           |       6 |

## Geometry (10 m pixels)

g/t: guide (10 m) / target (20 m) bands; gt_diff: difference of their global shifts;
rms_translation / rms_affine: block shifts left after one translation / an affine field.

|                   |   q10 |   q25 |   q50 |   q75 |   q90 |
|:------------------|------:|------:|------:|------:|------:|
| g_mag             | 0.125 | 0.172 | 0.278 | 0.554 | 0.788 |
| t_mag             | 0.126 | 0.197 | 0.331 | 0.601 | 0.907 |
| gt_diff           | 0.03  | 0.054 | 0.1   | 0.149 | 0.251 |
| rms_translation_g | 0.031 | 0.048 | 0.078 | 0.134 | 0.217 |
| rms_affine_g      | 0.024 | 0.037 | 0.063 | 0.104 | 0.169 |
| rms_translation_t | 0.076 | 0.101 | 0.174 | 0.268 | 0.464 |
| rms_affine_t      | 0.066 | 0.094 | 0.148 | 0.221 | 0.353 |

- tiles with |global shift| > 0.25 px: 58.5%
- tiles with |global shift| > 0.5 px: 28.2%
- tiles with |global shift| > 1.0 px: 5.5%

### H0 (one translation per tile leaves <= 0.25 px)

Pre-specified measure (rms of block shifts around their median):

- guide bands (10 m): available for 416 of 419 tiles; blocks per tile (median, same tiles): 16; median 0.078, q90 0.217; share <= 0.25: **91.6%** (site-bootstrap CI 87.0%-95.1%); missing counted as failing: 90.9% (CI 86.3%-94.4%)
- target bands (20 m): available for 415 of 419 tiles; blocks per tile (median, same tiles): 16; median 0.174, q90 0.464; share <= 0.25: **71.8%** (site-bootstrap CI 62.6%-79.6%); missing counted as failing: 71.1% (CI 61.9%-78.8%)

After the shifts the correction uses (guide bands: guide tile shift; target bands: joint target-band shift); rms of block shifts around 0:

- guide bands: available for 416 of 419 tiles; blocks per tile (median, same tiles): 16; median 0.082, q90 0.221; share <= 0.25: **92.3%** (site-bootstrap CI 87.7%-95.9%); missing counted as failing: 91.6% (CI 86.9%-95.3%)
- target bands (all six, raw; includes estimation noise): available for 416 of 419 tiles; blocks per tile (median, same tiles): 15; median 0.181, q90 0.457; share <= 0.25: **69.7%** (site-bootstrap CI 60.1%-77.8%); missing counted as failing: 69.2% (CI 59.5%-77.5%)
- target group red edge B5-B8A (raw): available for 416 of 419 tiles; median 0.237, q90 1.441; share <= 0.25: **52.2%** (site-bootstrap CI 40.3%-62.4%); missing counted as failing: 51.8% (CI 39.9%-62.0%)
- target group SWIR B11-B12 (raw): available for 408 of 419 tiles; median 0.256, q90 0.890; share <= 0.25: **48.8%** (site-bootstrap CI 36.5%-59.6%); missing counted as failing: 47.5% (CI 35.6%-58.3%)

H0' (post hoc): field common to the two target groups, unbiased squared estimate common2 = mean(a . b) thresholded at 0.25^2 (the clipped sqrt `common` is biased up by ~0.04 px at 14 blocks):

- H0' common target field (common2 <= 0.0625): available for 345 of 419 tiles; blocks per tile (median, same tiles): 15; median 0.003, q90 0.112; share <= 0.0625: **85.2%** (site-bootstrap CI 80.0%-89.6%); missing counted as failing: 70.2% (CI 63.8%-75.8%)
- field common to target and guide (common2 <= 0.0625): available for 345 of 419 tiles; blocks per tile (median, same tiles): 15; median 0.005, q90 0.052; share <= 0.0625: **91.3%** (site-bootstrap CI 87.3%-94.7%); missing counted as failing: 75.2% (CI 68.9%-80.9%)
- estimation noise of one target group (rms): available for 345 of 419 tiles; blocks per tile (median, same tiles): 15; median 0.258, q90 2.144; share <= 0.25: **48.7%** (site-bootstrap CI 36.9%-60.1%); missing counted as failing: 40.1% (CI 30.1%-49.1%)

Sensitivity: target bands aligned with the guide shift instead:

- target bands (raw): available for 416 of 419 tiles; median 0.206, q90 0.583; share <= 0.25: **62.5%** (site-bootstrap CI 51.9%-71.3%); missing counted as failing: 62.1% (CI 51.3%-70.9%)
- H0' measure: available for 345 of 419 tiles; median 0.010, q90 0.156; share <= 0.0625: **80.9%** (site-bootstrap CI 73.6%-86.7%); missing counted as failing: 66.6% (CI 58.8%-72.9%)

Sensitivity, excluding the top cloud-score quartile (45 tiles):

- guide bands: available for 372 of 374 tiles; median 0.080, q90 0.199; share <= 0.25: **94.9%** (site-bootstrap CI 92.2%-97.2%); missing counted as failing: 94.4% (CI 91.5%-96.7%)
- H0' measure: available for 310 of 374 tiles; median 0.003, q90 0.078; share <= 0.0625: **87.4%** (site-bootstrap CI 83.0%-91.1%); missing counted as failing: 72.5% (CI 66.8%-77.6%)

Geometry by |dt| (days): if residual fields were surface change they would grow with |dt| (R3-1):

| dt_bin   |   n |   rms_tileshift_g |   rms_tileshift_t |   common_t |   H0prime_share |   common_tg |
|:---------|----:|------------------:|------------------:|-----------:|----------------:|------------:|
| 0-3      |  40 |             0.049 |             0.103 |      0.044 |           0.968 |       0.039 |
| 4-7      |  36 |             0.061 |             0.101 |      0     |           1     |       0.04  |
| 8-15     |  45 |             0.081 |             0.2   |      0.01  |           0.8   |       0.083 |
| 16-40    | 298 |             0.092 |             0.205 |      0.068 |           0.826 |       0.078 |

Edge vs interior blocks after the tile shift (medians, 10 m px; reflect-padded tile edges, R3-3.NEW-3): rms_edge_g 0.085, rms_inner_g 0.064, rms_edge_t 0.185, rms_inner_t 0.150


Phase-correlation peaks:

|        |   q10 |   q25 |   q50 |   q75 |   q90 |
|:-------|------:|------:|------:|------:|------:|
| g_peak | 0.483 | 0.613 | 0.733 | 0.813 | 0.862 |
| t_peak | 0.469 | 0.592 | 0.716 | 0.803 | 0.87  |

### Where shifts vary (share of sum of squares; R1-1.M2)

|      |   site |   flight_in_site |   tile_in_flight |
|:-----|-------:|-----------------:|-----------------:|
| g_dy |  0.343 |            0.485 |            0.173 |
| g_dx |  0.354 |            0.424 |            0.222 |

## Radiometry per band, cloud-clean subset (374 tiles, cloud score <= q75; medians; reflectance units)

Statistics use only pixels whose whole 4-sigma low-pass footprint is valid. hp_gain is an *apparent* gain: S2 noise and aliasing in x lower it; not yet corrected (pre-E2 gate).

rmse_al: after global shift; rmse_lin_*: after shift + gain/offset fitted on the low-pass (ols, symmetric, or cross-fitted on the other tile half); hp_gain: std / OLS ratio of high-pass detail r vs x.

| band   |   rmse_raw |   rmse_al |   rmse_lin_ols |   rmse_lin_sym |   rmse_lin_cv |   a_ols |   a_sd |   r2_lowpass |   hp_gain_sd |   hp_gain_ols |   hp_corr_raw |   hp_corr |
|:-------|-----------:|----------:|---------------:|---------------:|--------------:|--------:|-------:|-------------:|-------------:|--------------:|--------------:|----------:|
| B02    |     0.0078 |    0.0073 |         0.0047 |         0.005  |        0.0053 |  0.674  | 0.7535 |       0.8113 |       0.5927 |        0.4843 |        0.7739 |    0.8311 |
| B03    |     0.0097 |    0.0092 |         0.0055 |         0.006  |        0.0062 |  0.6312 | 0.7281 |       0.8051 |       0.5629 |        0.4703 |        0.7601 |    0.8299 |
| B04    |     0.0129 |    0.012  |         0.0075 |         0.0079 |        0.0083 |  0.6818 | 0.7463 |       0.8572 |       0.6525 |        0.5563 |        0.8044 |    0.8544 |
| B08    |     0.0286 |    0.0279 |         0.0142 |         0.0152 |        0.0161 |  0.5139 | 0.6303 |       0.7302 |       0.5224 |        0.3971 |        0.7125 |    0.7857 |
| B05    |     0.0135 |    0.0133 |         0.0062 |         0.0065 |        0.0071 |  0.53   | 0.5972 |       0.8089 |       0.503  |        0.4111 |        0.7888 |    0.8245 |
| B06    |     0.0227 |    0.0223 |         0.009  |         0.0094 |        0.0103 |  0.4229 | 0.5213 |       0.7314 |       0.4312 |        0.3252 |        0.7214 |    0.7571 |
| B07    |     0.0257 |    0.0255 |         0.0106 |         0.0111 |        0.0122 |  0.4377 | 0.5345 |       0.7297 |       0.4464 |        0.3358 |        0.7164 |    0.7631 |
| B8A    |     0.0267 |    0.0264 |         0.011  |         0.0113 |        0.0125 |  0.439  | 0.5326 |       0.7451 |       0.4538 |        0.3483 |        0.7175 |    0.7627 |
| B11    |     0.0214 |    0.0212 |         0.0095 |         0.0101 |        0.0107 |  0.5042 | 0.5591 |       0.8203 |       0.5719 |        0.4895 |        0.8039 |    0.8412 |
| B12    |     0.0177 |    0.0172 |         0.0081 |         0.0086 |        0.0097 |  0.5422 | 0.6079 |       0.8456 |       0.6052 |        0.5153 |        0.8255 |    0.8607 |

RMSE reduction computed from the medians above (R1-3.7):

| band   |   align_pct (ratio of medians) |   lin_cv_pct (ratio of medians) |
|:-------|-------------------------------:|--------------------------------:|
| B02    |                            5.4 |                            28.4 |
| B03    |                            4.9 |                            33   |
| B04    |                            7.3 |                            30.8 |
| B08    |                            2.3 |                            42.3 |
| B05    |                            1.2 |                            46.5 |
| B06    |                            1.8 |                            53.9 |
| B07    |                            0.8 |                            52.2 |
| B8A    |                            1.1 |                            52.6 |
| B11    |                            1.2 |                            49.6 |
| B12    |                            2.7 |                            43.9 |

All reliable tiles incl. cloudy ones (419; sensitivity):

| band   |   rmse_lin_cv |   a_sd |   hp_gain_sd |   hp_corr |
|:-------|--------------:|-------:|-------------:|----------:|
| B02    |        0.0053 | 0.75   |       0.5817 |    0.8204 |
| B03    |        0.0062 | 0.7163 |       0.5517 |    0.8202 |
| B04    |        0.0083 | 0.7417 |       0.6463 |    0.8472 |
| B08    |        0.0164 | 0.6096 |       0.5063 |    0.77   |
| B05    |        0.0071 | 0.5918 |       0.4947 |    0.809  |
| B06    |        0.0105 | 0.4969 |       0.4169 |    0.7484 |
| B07    |        0.0124 | 0.5057 |       0.4294 |    0.7474 |
| B8A    |        0.0126 | 0.5082 |       0.4414 |    0.749  |
| B11    |        0.0107 | 0.5505 |       0.5662 |    0.8288 |
| B12    |        0.0095 | 0.5959 |       0.5937 |    0.8516 |

## Harmonisation check

A per-tile regression of S2 on the reference would match tile means exactly and shrink the reference's standard deviation to about R times that of S2 at the fitted scale.

| band   |   median R |   median sd ratio |   median hp_gain_sd |   corr(sd ratio, R) over tiles |   median |mean_r - mean_x| |
|:-------|-----------:|------------------:|--------------------:|-------------------------------:|---------------------------:|
| B02    |     0.8777 |            0.7104 |              0.5927 |                         0.7177 |                     0.0004 |
| B03    |     0.8652 |            0.684  |              0.5629 |                         0.7153 |                     0.0006 |
| B04    |     0.8957 |            0.7383 |              0.6525 |                         0.7497 |                     0.0008 |
| B08    |     0.8233 |            0.6079 |              0.5224 |                         0.7285 |                     0.0016 |
| B05    |     0.8746 |            0.5886 |              0.503  |                         0.7043 |                     0.0011 |
| B06    |     0.8259 |            0.5015 |              0.4312 |                         0.6486 |                     0.0017 |
| B07    |     0.832  |            0.5133 |              0.4464 |                         0.6868 |                     0.0021 |
| B8A    |     0.8364 |            0.5156 |              0.4538 |                         0.6639 |                     0.0022 |
| B11    |     0.8881 |            0.5702 |              0.5719 |                         0.6184 |                     0.002  |
| B12    |     0.9051 |            0.6143 |              0.6052 |                         0.6511 |                     0.0016 |

## Target bands: paired gain of each band's own shift over the joint target-band shift

d_rmse = rmse(tile shift) - rmse(own shift), d_hpcorr = hp_corr(own) - hp_corr(tile); positive = own shift better. Rows: tiles where the band's own peak >= min_peak.

| band   |   n |   own_shift_q50 |   own_shift_q90 |   d_rmse_median |   d_rmse_q90 |   share_rmse_gain_gt_0_001 |   d_hpcorr_median |   d_hpcorr_q90 |   share_hpcorr_gain_gt_0_01 |
|:-------|----:|----------------:|----------------:|----------------:|-------------:|---------------------------:|------------------:|---------------:|----------------------------:|
| B05    | 373 |          0.2953 |          0.8229 |               0 |       0      |                          0 |            0.0005 |         0.0069 |                      0.0509 |
| B06    | 372 |          0.3661 |          0.9484 |               0 |       0.0001 |                          0 |            0.0007 |         0.0073 |                      0.0753 |
| B07    | 372 |          0.3821 |          0.9815 |               0 |       0.0001 |                          0 |            0.001  |         0.0096 |                      0.0887 |
| B8A    | 372 |          0.396  |          0.9706 |               0 |       0.0001 |                          0 |            0.0014 |         0.0099 |                      0.0968 |
| B11    | 373 |          0.3578 |          0.9327 |               0 |       0.0001 |                          0 |            0.0027 |         0.01   |                      0.0992 |
| B12    | 373 |          0.3298 |          0.9101 |               0 |       0.0001 |                          0 |            0.0021 |         0.0102 |                      0.1046 |

By quartile of the guide-target shift difference gt_diff (10 m px):

| gt_diff         |   n |   d_rmse_median |   d_hpcorr_median |   share_hpcorr_gain_gt_0_01 |
|:----------------|----:|----------------:|------------------:|----------------------------:|
| (-0.001, 0.051] | 574 |               0 |            0.0009 |                      0.0645 |
| (0.051, 0.0943] | 552 |               0 |            0.001  |                      0.0634 |
| (0.0943, 0.141] | 558 |               0 |            0.0016 |                      0.0771 |
| (0.141, 1.023]  | 551 |               0 |            0.0022 |                      0.1397 |

## Cross-fitted residual RMSE by |dt| (days)

tiles per bin: 0-3: 40, 4-7: 36, 8-15: 45, 16-40: 298

| band   |    0-3 |    4-7 |   8-15 |   16-40 |
|:-------|-------:|-------:|-------:|--------:|
| B02    | 0.0052 | 0.0056 | 0.0044 |  0.0054 |
| B03    | 0.006  | 0.0068 | 0.005  |  0.0064 |
| B04    | 0.0065 | 0.0091 | 0.0068 |  0.009  |
| B08    | 0.014  | 0.0177 | 0.0168 |  0.0159 |
| B05    | 0.0058 | 0.0092 | 0.0056 |  0.0073 |
| B06    | 0.0086 | 0.0117 | 0.0105 |  0.01   |
| B07    | 0.0095 | 0.0142 | 0.0124 |  0.0119 |
| B8A    | 0.0098 | 0.0148 | 0.013  |  0.0123 |
| B11    | 0.0101 | 0.013  | 0.0093 |  0.0113 |
| B12    | 0.0079 | 0.0121 | 0.0072 |  0.0099 |

## Same, by land-cover superclass (medians; small classes are indicative only)

tiles per class: Rural: 212, Forest: 154, Water: 21, Developed: 7

| band   |   Developed |   Forest |   Rural |   Water |
|:-------|------------:|---------:|--------:|--------:|
| B02    |      0.0127 |   0.0043 |  0.0061 |  0.0044 |
| B03    |      0.0129 |   0.0053 |  0.007  |  0.0058 |
| B04    |      0.0141 |   0.0063 |  0.0109 |  0.0067 |
| B08    |      0.023  |   0.0164 |  0.0142 |  0.0318 |
| B05    |      0.0118 |   0.0055 |  0.0082 |  0.0075 |
| B06    |      0.0115 |   0.009  |  0.0103 |  0.0196 |
| B07    |      0.0144 |   0.0111 |  0.0119 |  0.0259 |
| B8A    |      0.0145 |   0.0117 |  0.0122 |  0.0271 |
| B11    |      0.0123 |   0.0095 |  0.012  |  0.0145 |
| B12    |      0.0161 |   0.0079 |  0.012  |  0.0094 |

## High-pass gain (sd ratio) of target bands by land cover

| band   |   Developed |   Forest |   Rural |   Water |
|:-------|------------:|---------:|--------:|--------:|
| B05    |       0.615 |    0.455 |   0.525 |   0.465 |
| B06    |       0.501 |    0.382 |   0.446 |   0.572 |
| B07    |       0.516 |    0.397 |   0.446 |   0.614 |
| B8A    |       0.556 |    0.415 |   0.456 |   0.626 |
| B11    |       0.574 |    0.559 |   0.579 |   0.74  |
| B12    |       0.652 |    0.594 |   0.605 |   0.666 |

## By tile cloud score (cloud_score_plus_cdf quartiles; SEN2NEON has no per-pixel cloud mask). Radiometry columns: medians over the six target bands and tiles.

| cloud_score    |   n |   g_peak |   rms_tileshift_g |   rmse_lin_cv |     a_sd |   hp_gain_sd |   hp_corr |
|:---------------|----:|---------:|------------------:|--------------:|---------:|-------------:|----------:|
| (0.881, 0.906] | 108 |   0.7689 |            0.0506 |        0.0075 |   0.5398 |       0.4744 |    0.8065 |
| (0.906, 0.914] | 129 |   0.7695 |            0.0814 |        0.0113 |   0.5803 |       0.5169 |    0.8221 |
| (0.914, 0.922] | 137 |   0.705  |            0.1057 |      nan      | nan      |     nan      |  nan      |
| (0.922, 0.937] |  45 |   0.497  |            0.1527 |      nan      | nan      |     nan      |  nan      |

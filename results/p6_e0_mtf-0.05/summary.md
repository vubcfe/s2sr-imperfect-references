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
| g_mag             | 0.125 | 0.172 | 0.278 | 0.552 | 0.794 |
| t_mag             | 0.128 | 0.19  | 0.326 | 0.6   | 0.907 |
| gt_diff           | 0.028 | 0.054 | 0.1   | 0.15  | 0.254 |
| rms_translation_g | 0.031 | 0.048 | 0.079 | 0.134 | 0.212 |
| rms_affine_g      | 0.025 | 0.037 | 0.064 | 0.103 | 0.167 |
| rms_translation_t | 0.075 | 0.102 | 0.175 | 0.263 | 0.445 |
| rms_affine_t      | 0.067 | 0.094 | 0.147 | 0.218 | 0.348 |

- tiles with |global shift| > 0.25 px: 58.5%
- tiles with |global shift| > 0.5 px: 28.2%
- tiles with |global shift| > 1.0 px: 5.5%

### H0 (one translation per tile leaves <= 0.25 px)

Pre-specified measure (rms of block shifts around their median):

- guide bands (10 m): available for 416 of 419 tiles; blocks per tile (median, same tiles): 16; median 0.079, q90 0.212; share <= 0.25: **91.8%** (site-bootstrap CI 87.3%-95.4%); missing counted as failing: 91.2% (CI 86.6%-94.7%)
- target bands (20 m): available for 414 of 419 tiles; blocks per tile (median, same tiles): 16; median 0.175, q90 0.445; share <= 0.25: **72.7%** (site-bootstrap CI 63.5%-80.3%); missing counted as failing: 71.8% (CI 62.6%-79.5%)

After the shifts the correction uses (guide bands: guide tile shift; target bands: joint target-band shift); rms of block shifts around 0:

- guide bands: available for 416 of 419 tiles; blocks per tile (median, same tiles): 16; median 0.083, q90 0.224; share <= 0.25: **91.6%** (site-bootstrap CI 86.8%-95.4%); missing counted as failing: 90.9% (CI 86.0%-94.7%)
- target bands (all six, raw; includes estimation noise): available for 416 of 419 tiles; blocks per tile (median, same tiles): 15; median 0.181, q90 0.459; share <= 0.25: **69.2%** (site-bootstrap CI 59.5%-77.2%); missing counted as failing: 68.7% (CI 58.9%-77.0%)
- target group red edge B5-B8A (raw): available for 416 of 419 tiles; median 0.236, q90 1.462; share <= 0.25: **51.7%** (site-bootstrap CI 39.2%-62.1%); missing counted as failing: 51.3% (CI 38.8%-61.7%)
- target group SWIR B11-B12 (raw): available for 408 of 419 tiles; median 0.257, q90 0.871; share <= 0.25: **48.3%** (site-bootstrap CI 35.8%-59.3%); missing counted as failing: 47.0% (CI 34.9%-58.0%)

H0' (post hoc): field common to the two target groups, unbiased squared estimate common2 = mean(a . b) thresholded at 0.25^2 (the clipped sqrt `common` is biased up by ~0.04 px at 14 blocks):

- H0' common target field (common2 <= 0.0625): available for 342 of 419 tiles; blocks per tile (median, same tiles): 14; median 0.003, q90 0.106; share <= 0.0625: **86.0%** (site-bootstrap CI 80.6%-90.6%); missing counted as failing: 70.2% (CI 63.5%-75.9%)
- field common to target and guide (common2 <= 0.0625): available for 342 of 419 tiles; blocks per tile (median, same tiles): 14; median 0.005, q90 0.047; share <= 0.0625: **92.4%** (site-bootstrap CI 88.3%-95.8%); missing counted as failing: 75.4% (CI 69.1%-81.3%)
- estimation noise of one target group (rms): available for 342 of 419 tiles; blocks per tile (median, same tiles): 14; median 0.253, q90 1.764; share <= 0.25: **49.4%** (site-bootstrap CI 37.3%-61.1%); missing counted as failing: 40.3% (CI 30.1%-49.5%)

Sensitivity: target bands aligned with the guide shift instead:

- target bands (raw): available for 416 of 419 tiles; median 0.208, q90 0.586; share <= 0.25: **62.5%** (site-bootstrap CI 52.1%-71.3%); missing counted as failing: 62.1% (CI 51.5%-71.0%)
- H0' measure: available for 341 of 419 tiles; median 0.010, q90 0.134; share <= 0.0625: **81.5%** (site-bootstrap CI 74.0%-87.6%); missing counted as failing: 66.3% (CI 58.1%-73.1%)

Sensitivity, excluding the top cloud-score quartile (45 tiles):

- guide bands: available for 372 of 374 tiles; median 0.080, q90 0.201; share <= 0.25: **94.1%** (site-bootstrap CI 91.3%-96.6%); missing counted as failing: 93.6% (CI 90.6%-96.2%)
- H0' measure: available for 309 of 374 tiles; median 0.003, q90 0.077; share <= 0.0625: **88.0%** (site-bootstrap CI 83.4%-91.9%); missing counted as failing: 72.7% (CI 66.9%-78.0%)

Geometry by |dt| (days): if residual fields were surface change they would grow with |dt| (R3-1):

| dt_bin   |   n |   rms_tileshift_g |   rms_tileshift_t |   common_t |   H0prime_share |   common_tg |
|:---------|----:|------------------:|------------------:|-----------:|----------------:|------------:|
| 0-3      |  40 |             0.051 |             0.101 |      0.049 |           0.968 |       0.041 |
| 4-7      |  36 |             0.063 |             0.102 |      0     |           1     |       0.038 |
| 8-15     |  45 |             0.084 |             0.195 |      0.005 |           0.829 |       0.081 |
| 16-40    | 298 |             0.094 |             0.205 |      0.07  |           0.832 |       0.076 |

Edge vs interior blocks after the tile shift (medians, 10 m px; reflect-padded tile edges, R3-3.NEW-3): rms_edge_g 0.087, rms_inner_g 0.065, rms_edge_t 0.187, rms_inner_t 0.149


Phase-correlation peaks:

|        |   q10 |   q25 |   q50 |   q75 |   q90 |
|:-------|------:|------:|------:|------:|------:|
| g_peak | 0.479 | 0.611 | 0.732 | 0.809 | 0.861 |
| t_peak | 0.47  | 0.593 | 0.718 | 0.801 | 0.868 |

### Where shifts vary (share of sum of squares; R1-1.M2)

|      |   site |   flight_in_site |   tile_in_flight |
|:-----|-------:|-----------------:|-----------------:|
| g_dy |  0.343 |            0.484 |            0.173 |
| g_dx |  0.354 |            0.424 |            0.222 |

## Radiometry per band, cloud-clean subset (374 tiles, cloud score <= q75; medians; reflectance units)

Statistics use only pixels whose whole 4-sigma low-pass footprint is valid. hp_gain is an *apparent* gain: S2 noise and aliasing in x lower it; not yet corrected (pre-E2 gate).

rmse_al: after global shift; rmse_lin_*: after shift + gain/offset fitted on the low-pass (ols, symmetric, or cross-fitted on the other tile half); hp_gain: std / OLS ratio of high-pass detail r vs x.

| band   |   rmse_raw |   rmse_al |   rmse_lin_ols |   rmse_lin_sym |   rmse_lin_cv |   a_ols |   a_sd |   r2_lowpass |   hp_gain_sd |   hp_gain_ols |   hp_corr_raw |   hp_corr |
|:-------|-----------:|----------:|---------------:|---------------:|--------------:|--------:|-------:|-------------:|-------------:|--------------:|--------------:|----------:|
| B02    |     0.0078 |    0.0075 |         0.0047 |         0.0051 |        0.0053 |  0.6738 | 0.7531 |       0.8113 |       0.5575 |        0.4487 |        0.7663 |    0.8201 |
| B03    |     0.0098 |    0.0094 |         0.0056 |         0.0059 |        0.0063 |  0.6311 | 0.7278 |       0.8051 |       0.5274 |        0.4364 |        0.7485 |    0.8157 |
| B04    |     0.0131 |    0.0122 |         0.0076 |         0.008  |        0.0083 |  0.6815 | 0.746  |       0.8573 |       0.6139 |        0.5178 |        0.8003 |    0.8466 |
| B08    |     0.0287 |    0.0282 |         0.014  |         0.0152 |        0.016  |  0.5136 | 0.63   |       0.7302 |       0.4867 |        0.3611 |        0.7113 |    0.7777 |
| B05    |     0.0137 |    0.0136 |         0.0062 |         0.0065 |        0.007  |  0.5297 | 0.5962 |       0.8088 |       0.449  |        0.3667 |        0.7886 |    0.8199 |
| B06    |     0.0226 |    0.0224 |         0.0089 |         0.0092 |        0.0101 |  0.4227 | 0.5207 |       0.7314 |       0.383  |        0.2884 |        0.7214 |    0.7544 |
| B07    |     0.0257 |    0.0256 |         0.0103 |         0.0109 |        0.012  |  0.4371 | 0.5338 |       0.7295 |       0.398  |        0.3014 |        0.718  |    0.7582 |
| B8A    |     0.0268 |    0.0266 |         0.0107 |         0.0111 |        0.0122 |  0.4389 | 0.5321 |       0.7454 |       0.4068 |        0.3079 |        0.7208 |    0.7601 |
| B11    |     0.0216 |    0.0215 |         0.0092 |         0.01   |        0.0105 |  0.5038 | 0.5577 |       0.8203 |       0.4915 |        0.4246 |        0.8073 |    0.8425 |
| B12    |     0.018  |    0.0177 |         0.008  |         0.0085 |        0.0094 |  0.5413 | 0.6074 |       0.8461 |       0.5205 |        0.4439 |        0.8286 |    0.859  |

RMSE reduction computed from the medians above (R1-3.7):

| band   |   align_pct (ratio of medians) |   lin_cv_pct (ratio of medians) |
|:-------|-------------------------------:|--------------------------------:|
| B02    |                            3.4 |                            28.7 |
| B03    |                            4.2 |                            33.5 |
| B04    |                            6.2 |                            32.3 |
| B08    |                            1.8 |                            43.1 |
| B05    |                            0.8 |                            48.1 |
| B06    |                            0.6 |                            55.2 |
| B07    |                            0.3 |                            53.2 |
| B8A    |                            0.9 |                            54.1 |
| B11    |                            0.6 |                            51.2 |
| B12    |                            1.8 |                            46.7 |

All reliable tiles incl. cloudy ones (419; sensitivity):

| band   |   rmse_lin_cv |   a_sd |   hp_gain_sd |   hp_corr |
|:-------|--------------:|-------:|-------------:|----------:|
| B02    |        0.0053 | 0.7497 |       0.5475 |    0.8099 |
| B03    |        0.0062 | 0.7161 |       0.5163 |    0.8066 |
| B04    |        0.0083 | 0.7413 |       0.6081 |    0.8375 |
| B08    |        0.0163 | 0.6094 |       0.4709 |    0.7618 |
| B05    |        0.007  | 0.5909 |       0.4416 |    0.8071 |
| B06    |        0.0103 | 0.4964 |       0.3713 |    0.7469 |
| B07    |        0.0122 | 0.5049 |       0.3827 |    0.7456 |
| B8A    |        0.0125 | 0.5076 |       0.3933 |    0.7515 |
| B11    |        0.0105 | 0.5492 |       0.4855 |    0.8288 |
| B12    |        0.0093 | 0.5948 |       0.5125 |    0.8488 |

## Harmonisation check

A per-tile regression of S2 on the reference would match tile means exactly and shrink the reference's standard deviation to about R times that of S2 at the fitted scale.

| band   |   median R |   median sd ratio |   median hp_gain_sd |   corr(sd ratio, R) over tiles |   median |mean_r - mean_x| |
|:-------|-----------:|------------------:|--------------------:|-------------------------------:|---------------------------:|
| B02    |     0.8745 |            0.6999 |              0.5575 |                         0.7158 |                     0.0004 |
| B03    |     0.86   |            0.6727 |              0.5274 |                         0.7152 |                     0.0006 |
| B04    |     0.8965 |            0.7295 |              0.6139 |                         0.7478 |                     0.0008 |
| B08    |     0.8217 |            0.5963 |              0.4867 |                         0.7271 |                     0.0016 |
| B05    |     0.8742 |            0.5758 |              0.449  |                         0.6977 |                     0.0011 |
| B06    |     0.8224 |            0.4924 |              0.383  |                         0.6424 |                     0.0017 |
| B07    |     0.8302 |            0.5036 |              0.398  |                         0.6808 |                     0.0021 |
| B8A    |     0.8337 |            0.5084 |              0.4068 |                         0.6578 |                     0.0022 |
| B11    |     0.89   |            0.5553 |              0.4915 |                         0.6087 |                     0.002  |
| B12    |     0.9056 |            0.5973 |              0.5205 |                         0.6408 |                     0.0016 |

## Target bands: paired gain of each band's own shift over the joint target-band shift

d_rmse = rmse(tile shift) - rmse(own shift), d_hpcorr = hp_corr(own) - hp_corr(tile); positive = own shift better. Rows: tiles where the band's own peak >= min_peak.

| band   |   n |   own_shift_q50 |   own_shift_q90 |   d_rmse_median |   d_rmse_q90 |   share_rmse_gain_gt_0_001 |   d_hpcorr_median |   d_hpcorr_q90 |   share_hpcorr_gain_gt_0_01 |
|:-------|----:|----------------:|----------------:|----------------:|-------------:|---------------------------:|------------------:|---------------:|----------------------------:|
| B05    | 373 |          0.2953 |          0.8229 |               0 |       0      |                          0 |            0.0007 |         0.0063 |                      0.0456 |
| B06    | 372 |          0.3674 |          0.9469 |               0 |       0.0001 |                          0 |            0.0008 |         0.007  |                      0.0538 |
| B07    | 372 |          0.3821 |          0.9817 |               0 |       0.0001 |                          0 |            0.001  |         0.0081 |                      0.0726 |
| B8A    | 372 |          0.396  |          0.9706 |               0 |       0.0001 |                          0 |            0.0012 |         0.0087 |                      0.0753 |
| B11    | 373 |          0.3561 |          0.9298 |               0 |       0.0001 |                          0 |            0.0024 |         0.0089 |                      0.0724 |
| B12    | 373 |          0.328  |          0.9101 |               0 |       0.0001 |                          0 |            0.002  |         0.0103 |                      0.1046 |

By quartile of the guide-target shift difference gt_diff (10 m px):

| gt_diff         |   n |   d_rmse_median |   d_hpcorr_median |   share_hpcorr_gain_gt_0_01 |
|:----------------|----:|----------------:|------------------:|----------------------------:|
| (-0.001, 0.051] | 574 |               0 |            0.0007 |                      0.054  |
| (0.051, 0.0922] | 546 |               0 |            0.0009 |                      0.0458 |
| (0.0922, 0.143] | 588 |               0 |            0.0014 |                      0.0595 |
| (0.143, 1.113]  | 527 |               0 |            0.002  |                      0.1271 |

## Cross-fitted residual RMSE by |dt| (days)

tiles per bin: 0-3: 40, 4-7: 36, 8-15: 45, 16-40: 298

| band   |    0-3 |    4-7 |   8-15 |   16-40 |
|:-------|-------:|-------:|-------:|--------:|
| B02    | 0.0053 | 0.0059 | 0.0045 |  0.0055 |
| B03    | 0.006  | 0.0069 | 0.005  |  0.0064 |
| B04    | 0.0069 | 0.0091 | 0.0068 |  0.0089 |
| B08    | 0.0139 | 0.0177 | 0.0167 |  0.0158 |
| B05    | 0.0058 | 0.0091 | 0.0056 |  0.0072 |
| B06    | 0.0082 | 0.0117 | 0.0104 |  0.0099 |
| B07    | 0.0093 | 0.0138 | 0.0123 |  0.0119 |
| B8A    | 0.0098 | 0.0146 | 0.0128 |  0.0121 |
| B11    | 0.0099 | 0.013  | 0.0091 |  0.0111 |
| B12    | 0.0079 | 0.0119 | 0.0071 |  0.0097 |

## Same, by land-cover superclass (medians; small classes are indicative only)

tiles per class: Rural: 212, Forest: 154, Water: 21, Developed: 7

| band   |   Developed |   Forest |   Rural |   Water |
|:-------|------------:|---------:|--------:|--------:|
| B02    |      0.0133 |   0.0044 |  0.0061 |  0.0043 |
| B03    |      0.0131 |   0.0053 |  0.0071 |  0.0058 |
| B04    |      0.0143 |   0.0062 |  0.0108 |  0.0067 |
| B08    |      0.0233 |   0.0165 |  0.0143 |  0.0311 |
| B05    |      0.0119 |   0.0054 |  0.0081 |  0.0075 |
| B06    |      0.0113 |   0.009  |  0.0101 |  0.0193 |
| B07    |      0.0141 |   0.011  |  0.0118 |  0.0254 |
| B8A    |      0.0141 |   0.0115 |  0.012  |  0.027  |
| B11    |      0.0119 |   0.0092 |  0.0117 |  0.0145 |
| B12    |      0.0162 |   0.0077 |  0.0119 |  0.0093 |

## High-pass gain (sd ratio) of target bands by land cover

| band   |   Developed |   Forest |   Rural |   Water |
|:-------|------------:|---------:|--------:|--------:|
| B05    |       0.555 |    0.411 |   0.47  |   0.413 |
| B06    |       0.447 |    0.34  |   0.403 |   0.502 |
| B07    |       0.452 |    0.353 |   0.401 |   0.54  |
| B8A    |       0.487 |    0.366 |   0.409 |   0.551 |
| B11    |       0.495 |    0.476 |   0.493 |   0.642 |
| B12    |       0.558 |    0.514 |   0.52  |   0.58  |

## By tile cloud score (cloud_score_plus_cdf quartiles; SEN2NEON has no per-pixel cloud mask). Radiometry columns: medians over the six target bands and tiles.

| cloud_score    |   n |   g_peak |   rms_tileshift_g |   rmse_lin_cv |     a_sd |   hp_gain_sd |   hp_corr |
|:---------------|----:|---------:|------------------:|--------------:|---------:|-------------:|----------:|
| (0.881, 0.906] | 108 |   0.7653 |            0.0512 |        0.0075 |   0.539  |       0.4188 |    0.8048 |
| (0.906, 0.914] | 129 |   0.7655 |            0.0813 |        0.0111 |   0.5791 |       0.4539 |    0.8199 |
| (0.914, 0.922] | 137 |   0.7057 |            0.1059 |      nan      | nan      |     nan      |  nan      |
| (0.922, 0.937] |  45 |   0.4933 |            0.1527 |      nan      | nan      |     nan      |  nan      |

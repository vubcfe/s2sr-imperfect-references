# Harmonisation-regression check at 2.5 m (the SEN2NEON fitting grid)

Source: results/p6_e0/harmonisation_regression_check.csv (one row per tile and band). Per tile: R = correlation of the Sentinel-2 band, bilinearly interpolated to 2.5 m, with the released reference; sd_ratio = sd(reference)/sd(interpolated Sentinel-2); gap = sd_ratio - R. An OLS fit of Sentinel-2 on the airborne image predicts sd_ratio = R and equal tile means.
This differs from the 'Harmonisation check' table of summary.md, which compares the reference with native Sentinel-2 at 20 m / 10 m (a different scale, where the prediction sd_ratio = R does not apply).

| band   |   n_tiles |   median_R |   median_sd_ratio |   median_gap |   median_abs_gap |   corr_over_tiles |   median_abs_mean_diff |
|:-------|----------:|-----------:|------------------:|-------------:|-----------------:|------------------:|-----------------------:|
| B02    |      2100 |     0.8231 |            0.8185 |      -0.0179 |           0.0287 |            0.7734 |                 0.0001 |
| B04    |      2100 |     0.8402 |            0.8276 |      -0.0223 |           0.0311 |            0.7323 |                 0.0001 |
| B08    |      2100 |     0.7185 |            0.7077 |      -0.0241 |           0.0402 |            0.7911 |                 0.0001 |
| B05    |      2100 |     0.76   |            0.7286 |      -0.0323 |           0.0453 |            0.7496 |                 0.0001 |
| B06    |      2100 |     0.6741 |            0.6495 |      -0.0355 |           0.0504 |            0.7783 |                 0.0001 |
| B07    |      2100 |     0.6842 |            0.6581 |      -0.0341 |           0.0481 |            0.8068 |                 0.0001 |
| B8A    |      2100 |     0.6851 |            0.6553 |      -0.0354 |           0.05   |            0.7998 |                 0.0001 |
| B11    |      2100 |     0.7525 |            0.7113 |      -0.0423 |           0.056  |            0.7494 |                 0.0001 |
| B12    |      2100 |     0.7889 |            0.7501 |      -0.0388 |           0.0499 |            0.7457 |                 0.0001 |

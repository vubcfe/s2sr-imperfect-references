# Licence of the data and results

The results (`results/`, `runs/`), the figures (`paper/figures/`) and the derived data archived on Zenodo
(registered and corrected references, tile shifts, compression factors, site split, model weights) are released
under the **Creative Commons Attribution 4.0 International licence (CC-BY-4.0)**,
https://creativecommons.org/licenses/by/4.0/.

They are derived from **SEN2NEON**, which is distributed under CC-BY-4.0 by its authors:

> Donike, S., Aybar, C., Contreras, J., Gómez-Chova, L. (2026). SEN2NEON: Enabling quantitative benchmarking of
> Sentinel-2 superresolution for all multispectral bands. IEEE Geoscience and Remote Sensing Letters 23, 6013905.
> https://doi.org/10.1109/LGRS.2026.3703947 — data: https://huggingface.co/datasets/isp-uv-es/SEN2NEON

SEN2NEON contains modified Copernicus Sentinel-2 data and airborne imagery of the NEON Airborne Observation Platform
(AVIRIS-NG). When you use the derived data, cite this repository/paper and SEN2NEON, and keep the attributions of
Copernicus Sentinel data and NEON required by their providers.

The changes made to SEN2NEON are described in the paper and implemented in `scripts/p6_build_domain_r.py`
(registration by the per-tile shift of `results/p6_e0/tiles.csv` and one compression factor per tile and band).

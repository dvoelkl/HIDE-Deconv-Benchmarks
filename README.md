# Additional Benchmarks for HIDE-Deconv

This repository contains reproducible scripts for additional benchmarks for the HIDE-Deconv Bulk RNA-seq Deconvolution package:

**[HIDE-Deconv](https://github.com/dvoelkl/HIDE-deconv)**

## Benchmarks on simulated bulks (HIDE-Deconv v0.1.0)
For information on how the benchmarks on simulated lung adenocarcinoma bulks where performed, we refer to the preprint: 

HIDE-Deconv: A hierarchical deconvolution framework for multiscale characterization of cellular remodeling

Dennis Voelkl, Sarah Bolz, Austin Rayford, Thomas Stevenson, Thomas Sterr, Malte Mensching-Buhr, Nicole Seifert, Julia Arp, Cornelia Schuster, Jana Tausche, Laurenz Engel, Helena U. Zacharias, Michael Altenbuchinger, Franziska Görter bioRxiv 2026.08.24.746754; doi: https://doi.org/10.64898/2026.08.24.746754

## Benchmarks on real world data: Morandini (HIDE-Deconv v0.4.0)
As an additional evaluation on real world data after posting the preprint, HIDE-deconv was benchmarked on the Morandini bulk RNA-seq dataset using the HaoSub single-cell reference, following the evaluation procedure of the `deconvBench` framework from Dietrich et al. (2026).

| Method          | Pearson R   |    RMSE   |
| --------------- | ----------: | --------: |
| AutoGeneS       |       0.819 |     0.182 |
| BayesPrism      |       0.254 |     0.295 |
| Bisque          |       0.872 |     0.122 |
| CIBERSORTx      |       0.648 |     0.193 |
| DWLS            |       0.909 |     0.103 |
| MuSiC           |       0.319 |     0.275 |
| Scaden          |       0.795 |     0.150 |
| SCDC            |       0.468 |     0.239 |
| **HIDE-deconv (1)** |   **0.934** | **0.094** |
| **HIDE-deconv (2)** |   **0.933** | **0.094** |

[HIDE-deconv (1)](./morandini_benchmark_hidedeconv.py) was trained and evaluated at the resolution of the cell types defined in the Hao single-cell reference, without specifying a hierarchy.

[HIDE-deconv (2)](./morandini_benchmark_hidedeconv_wHierarchy.py) was trained and evaluated using a three-level cell-type hierarchy, analogous to the hierarchical aggregation used for the FACS benchmarking.

Both HIDE-Deconv runs were evaluated independently. Results for the other methods are taken from [Dietrich et al. (2026)](https://doi.org/10.1186/s13059-026-03955-w)

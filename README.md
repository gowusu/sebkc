# Surface Energy Balance and Crop Coefficient Estimation with R

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/sebkc)](https://CRAN.R-project.org/package=sebkc)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/sebkc)](https://CRAN.R-project.org/package=sebkc)
[![DOI](https://img.shields.io/badge/DOI-10.32614%2FCRAN.package.sebkc-blue)](https://doi.org/10.32614/CRAN.package.sebkc)
<!-- badges: end -->

This R package computes and integrates surface energy balance into the FAO56 crop water balance model.

Specifically, the package can perform the following functions:

+   Single crop coefficient modelling
+   Dual crop coefficient modelling
+   The integration of thermal-based Evaporative Fractions in a water balance model
+   The computation of Two-Source Surface Energy Balance (TSEB) models such as TSEB-PT,
    TSEB-PM and TSEB-Parallel
+   The computation of One-Source Surface Energy Balance (OSEB) models such as SEBAL, METRIC, SEBI, SSEB and SEBS

## Installation

Install the released version from CRAN (recommended — dependencies are handled automatically):

```r
install.packages("sebkc")
```

Or the development version from GitHub:

```r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("gowusu/sebkc")
```

## Citation

If you use **sebkc** in your work, please cite it. You can always get the current
citation in R with:

```r
citation("sebkc")
```

Owusu, G. (2026). *sebkc: Surface Energy Balance and Crop Coefficient Evapotranspiration
Estimation.* R package version 1.0-6. https://doi.org/10.32614/CRAN.package.sebkc

BibTeX:

```bibtex
@Manual{sebkc,
  title  = {sebkc: Surface Energy Balance and Crop Coefficient Evapotranspiration Estimation},
  author = {George Owusu},
  year   = {2026},
  note   = {R package version 1.0-6},
  doi    = {10.32614/CRAN.package.sebkc},
  url    = {https://doi.org/10.32614/CRAN.package.sebkc},
}
```

## Links

+   CRAN: https://CRAN.R-project.org/package=sebkc
+   DOI: https://doi.org/10.32614/CRAN.package.sebkc

**Author**: George Owusu <owusugeorge@ug.edu.gh>. This is part of the author's PhD program.

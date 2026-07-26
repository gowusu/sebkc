# Surface Energy Balance and Crop Coefficient Estimation with R

This R package computes and integrates surface energy balance into the
FAO56 crop water balance model.

Specifically, the package can perform the following functions:

- Single crop coefficient modelling
- Dual crop coefficient modelling
- The integration of thermal-based Evaporative Fractions in a water
  balance model
- The computation of Two-Source Surface Energy Balance (TSEB) models
  such as TSEB-PT, TSEB-PM and TSEB-Parallel
- The computation of One-Source Surface Energy Balance (OSEB) models
  such as SEBAL, METRIC, SEBI, SSEB and SEBS

## Installation

Install the released version from CRAN (recommended — dependencies are
handled automatically):

``` r
install.packages("sebkc")
```

Or the development version from GitHub:

``` r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("gowusu/sebkc")
```

## Documentation website

Full documentation — both manuals plus a reference for every function —
is online at **<https://gowusu.github.io/sebkc/>**.

## Manuals

Two manuals ship with the package — a plain-language guide and a full
technical reference — with every number and figure produced by running
the package.

**Start here if you just want crop water numbers —** *How Much Water
Does My Crop Need?* A practical, jargon-free guide that states a real
problem (how much to irrigate crops on farms near Kumasi) and walks
through the answer for maize, tomato, cabbage, onion and pepper across
the dry and main seasons, using the built-in Kumasi weather (kept up to
date automatically):

- **[Read it
  online](https://gowusu.github.io/sebkc/articles/crop_water_needs.html)**
  — on the documentation website
- [Source on GitHub
  (Markdown)](https://github.com/gowusu/sebkc/blob/master/manual/crop_water_needs.md)
  — renders in the repo
- `manual/crop_water_needs.docx` — Word version (Save As → PDF for a
  printable copy)
- `manual/reproduce_water_needs.R` — regenerates every table and figure

**The full technical reference —** *Evapotranspiration from Weather and
Satellites.* Reference evapotranspiration and the FAO-56
crop-coefficient water balance, plus the satellite
surface-energy-balance models (SEBAL, METRIC, SEBS, SEBI, SSEB, S-SEBI
and the two-source TSEB), on the bundled Kumasi Landsat data:

- **[Read it
  online](https://gowusu.github.io/sebkc/articles/sebkc_manual.html)** —
  on the documentation website
- [Source on GitHub
  (Markdown)](https://github.com/gowusu/sebkc/blob/master/manual/sebkc_manual.md)
  — renders in the repo
- `manual/sebkc_manual.docx` — Word version (Save As → PDF for a
  printable copy)
- `manual/reproduce_manual.R` — regenerates every table and figure

## Citation

If you use **sebkc** in your work, please cite it. You can always get
the current citation in R with:

``` r
citation("sebkc")
```

Owusu, G. (2026). *sebkc: Surface Energy Balance and Crop Coefficient
Evapotranspiration Estimation.* R package version 1.0-6.
<https://doi.org/10.32614/CRAN.package.sebkc>

BibTeX:

``` bibtex
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

- CRAN: <https://CRAN.R-project.org/package=sebkc>
- DOI: <https://doi.org/10.32614/CRAN.package.sebkc>

**Author**: George Owusu <owusugeorge@ug.edu.gh>. This is part of the
author’s PhD program.

Surface Energy Balance and Crop Coefficient Estimation with R
=============================================================

This R package computes and integrates surface energy balance into FAO56 crop water balance model

Specifically, the package can perform the following functions:

-   Single crop cofficient modelling
-   Dual crop cofficient modelling
-   the integration of thermal-based Evaporative Fractions in a water balance model
-   The computation of Two Source Surface Energy Balance (TSEB) model such as such as TSEB-PT, TSEB-PM and TSEB-Parallel
-   The computation of One Source Surface Energy Balance (OSEB) models such as SEBAL, METRIC, SEBI, SSEB and SEBS

Installation
------------

You can install the latest development version from github with,
<pre><code> if (packageVersion("devtools") < 1.6) {
  install.packages("devtools")
  }
  if(!require(sp)){
  install.packages("sp")
}
if(!require(raster)){
  install.packages("raster")
}
if(!require(gstat)){
  install.packages("gstat")
}
  devtools::install_github("gowusu/sebkc")
</code></pre>
**Author**: George Owusu <owusugeorge@ug.edu.gh> This is part of the author's PhD program.

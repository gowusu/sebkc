## Submission summary

This is a resubmission / first CRAN submission of sebkc 1.0-2.

sebkc provides surface energy balance and crop coefficient evapotranspiration
models (SEBAL, METRIC, SEBI, SSEB, SEBS, TSEB) integrated with the FAO56 water
balance model.

## Test environments

<!-- Fill these in with YOUR actual results before submitting -->
* local Windows install, R 4.x.x
* win-builder (devel and release)   <-- run devtools::check_win_devel() / check_win_release()
* macOS / Linux via R-hub           <-- run rhub::rhub_check()

## R CMD check results

<!-- Paste the real summary line here. Target: 0 errors | 0 warnings | 0 notes -->
0 errors | 0 warnings | 0 notes

## Notes for CRAN reviewers

* Examples that require interactive input (drawing an area of interest on a
  map via drawExtent()/drawPoly()), internet access (weather station / NASA
  downloads), or long run times are wrapped in \dontrun{}.
* The package writes output only to a user-supplied folder argument; no files
  are written to the user's home, working directory, or any location outside
  tempdir() during checks/examples.

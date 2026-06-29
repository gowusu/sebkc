## Submission summary

This is a new submission of sebkc 1.0-2.

sebkc provides surface energy balance and crop coefficient evapotranspiration
models (SEBAL, METRIC, SEBI, SSEB, SEBS, TSEB) integrated with the FAO56 water
balance model.

## Test environments

* local: Windows 10, R 4.6.0 -- 0 errors | 0 warnings | 0 notes
* win-builder R-devel and R-release

## R CMD check results

0 errors | 0 warnings | 2 notes

### NOTE 1: CRAN incoming feasibility

* New submission. (Expected.)
* "Possibly misspelled words in DESCRIPTION": these are not misspellings.
  They are model acronyms (SEBAL, METRIC, OSEB, SEBI, SEBS, SSEB, TSEB),
  the surname of a cited author (Bastiaanssen), the abbreviation "et al.",
  the organisation name FAO, and standard domain terms (evapotranspiration,
  evaporative).

### NOTE 2 (HTML manual math), now fixed

* The two flagged expressions in tseb were changed from \eqn{} to \code{}
  because they are R expressions (containing `$`), not mathematics.

## Other notes for reviewers

* Examples that require interactive input (drawing an area of interest on a
  map via drawExtent()/drawPoly()), internet access (weather-station / NASA
  downloads) or long run times are wrapped in \dontrun{}.
* The package writes output only to a user-supplied folder argument; nothing
  is written outside tempdir() during checks/examples.
* installed size is ~5.5 Mb because inst/extdata contains a small Landsat
  scene used by the runnable examples and required to demonstrate the models.

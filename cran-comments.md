## Submission summary

Resubmission of sebkc 1.0-3, addressing the review by Konstanze Lauseker.

## Changes made in response to the review

* Explained all acronyms (SEBAL, METRIC, SEBI, SSEB, SEBS, TSEB-PT/PM,
  OSEB, FAO-56) in the Description.
* Replaced `T`/`F` with `TRUE`/`FALSE` throughout.
* Replaced console `print()` output with `message()`.
* Removed `options(warn=-1)`.
* Save and restore `par()` via `on.exit()` in all functions that change it.
* Save and restore the user's RNG state around the internal `set.seed()`.
* Fixed unexecutable example code (SEBS, weather).
* Fixed a real length-1 bug in `biomass()` (`if(is.na(Ndays))` ->
  `if(any(is.na(Ndays)))`).

## On \dontrun vs \donttest

The pure-numeric FAO-56 reference examples (`ETo()`, `ETohr()`) are now in
`\donttest{}` and run during checks.

The remaining examples are kept in `\dontrun{}` because they genuinely
cannot run within check constraints:

* `weather()` downloads live weather-station / reanalysis data over the
  internet.
* The surface-energy-balance and crop-coefficient models (SEBAL, METRIC,
  SEBI, SSEB, SEBS, TSEB, kc, biomass, landsat578, hotTs, coldTs, ...)
  each process a full Landsat scene (and, for some, perform spatial
  interpolation), taking minutes per example and requiring the bundled
  ~5 MB scene. These are impractical to execute during routine checks.

## Test environments

* local: Windows 10, R 4.6.0 -- 0 errors | 0 warnings | 0 notes
* win-builder R-devel and R-release (1.0-2 baseline): 1 NOTE (New submission)

## R CMD check results

0 errors | 0 warnings | 1 note (New submission; installed size ~5.5 Mb due
to the bundled Landsat scene used by examples).

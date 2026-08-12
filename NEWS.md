# sebkc — NEWS

## Development version (branch fix-sebs-ef)

### Bug fixes
* **SEBAL/METRIC** (`sebal()`) no longer errors when hot/cold anchor pixels are
  supplied as coordinate pairs (or as `hotTs()`/`coldTs()` objects). The keyword
  test (`xyhot == "auto"` / `"full"`) was evaluated on the whole argument, which
  under R >= 4.3 raises *"'length = 2' in coercion to 'logical(1)'"* for any
  length-2 numeric anchor. The check now tests `is.numeric()` first and matches
  the keyword on a single element. This had been blocking every manual-anchor run
  (and therefore the full validation pipeline) on current R.
* **SEBI** (`sebi()`) now derives the dry and wet evaporative limits **per pixel**
  (location-dependent, after Menenti & Choudhury 1993) instead of from the single
  hot/cold anchor pixels, and corrects a sign error in the available-energy term
  (`Rn - G`, previously `Rn + G`). This removes the spurious dependence on the
  chosen spatial domain.

### Changes
* **SEBS** (`sebs()`) has been reformulated to compute its **own** sensible heat
  from the actual surface-to-air temperature difference and a roughness length for
  heat derived from Su (2002)'s kB^-1 excess-resistance parameterisation, with a
  Monin-Obukhov stability iteration. The previous version reused SEBAL's
  anchor-calibrated H and rescaled it by the cold-pixel latent heat, which (a) made
  SEBS depend on the chosen anchor pixels and (b) was algebraically redundant with
  SEBAL. The new implementation is domain-independent by construction.
  * **EXPERIMENTAL — validated and found wanting.** The reformulated H was
    validated against AmeriFlux eddy-covariance sensible heat at the Mead maize
    sites (US-Ne1, US-Ne2) for three cloud-free 2001 Landsat 7 overpasses
    (DOY 161/225/273, WRS-2 028/031). The native H does **not** reproduce the
    measured flux: it under-predicts H over dry/sparse canopy and its seasonal
    ordering is uncorrelated with observations (rank r ~ 0; RMSE ~150 W/m2). The
    cause is structural — the radiometric surface temperature diverges from the
    aerodynamic temperature over partial and senescent canopies (e.g. late season,
    measured H ~200 W/m2 while Ts < Ta), which a one-source model cannot represent —
    so it is **not** a scalar-kB^-1 miscalibration that could be tuned out. `sebs()`
    therefore remains experimental and continues to emit a warning. Use
    `sebal()`/METRIC (anchor-calibrated dT) or `tseb()` (two-source) for
    quantitative work.
  * `sebs()` gains a `Ta` argument for the actual overpass air temperature. The
    previous code used the daily `Tmax` as the air temperature, which biases the
    surface-to-air gradient (and hence H) low; supplying `Ta` is now recommended,
    and a message is emitted when it is omitted.

### Housekeeping
* Removed dead, machine-specific absolute-path comments from `R/ETo.R` and
  `R/ETohr.R` (portability).

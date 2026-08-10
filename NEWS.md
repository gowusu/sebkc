# sebkc — NEWS

## Development version (branch fix-sebs-ef)

### Bug fixes
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
  * **EXPERIMENTAL:** the kB^-1 excess resistance has not yet been calibrated or
    validated against measured sensible-heat flux. On current test data the absolute
    fluxes are not reliable (EF tends to saturate). `sebs()` now emits a warning to
    this effect. Use `sebal()`/METRIC for quantitative work until SEBS is validated.

### Housekeeping
* Removed dead, machine-specific absolute-path comments from `R/ETo.R` and
  `R/ETohr.R` (portability).

# CRAN readiness for `sebkc` — status & action plan

_Last updated: 2026-06-29_

CRAN acceptance is **defined** by `R CMD check --as-cran` returning
**0 errors, 0 warnings, 0 notes**, run iteratively until clean. That check
must be run in a real R + toolchain environment — it could not be run here.
This file lists what was already done and exactly what remains.

---

## ✅ Done (committed)

| Item | What changed |
|------|--------------|
| Runtime crash | Fixed `Error in -i` in `hotTs()`/`coldTs()` (NA-free rasters). |
| DESCRIPTION | `Authors@R`, Title Case title, reworded Description with `<doi:>` refs, `Encoding: UTF-8`, `URL`/`BugReports`, dropped `LazyData` (no `data/` dir), bumped `RoxygenNote`, declared base-package Imports. |
| `.Rbuildignore` | Ignore `README.html`, `README.RMD`, `cran-comments.md`, etc. |
| `cran-comments.md` | Template to fill with real check results. |

---

## ⛔ Must do before submitting (requires R on your machine)

Install the toolchain first:
```r
install.packages(c("devtools","roxygen2","rhub","spelling"))
# Windows: also install Rtools matching your R version (rtools is needed to
# BUILD/CHECK even pure-R packages on Windows for the source tarball step).
```

### 1. Regenerate documentation & NAMESPACE
The old code uses `@import raster` / `@import sp`. To remove the inevitable
"no visible global function" NOTEs you must add `@importFrom` tags for base
functions (e.g. `stats::kmeans`, `stats::quantile`, `graphics::hist`,
`graphics::par`, `graphics::abline`, `utils::read.csv`, `utils::read.table`,
`grDevices::dev.new`) and then:
```r
devtools::document()   # regenerates NAMESPACE + man/ from roxygen
```

### 2. Run the check loop until clean
```r
devtools::check()                # local --as-cran check
# fix every ERROR, WARNING, NOTE, then repeat
```

### 3. Cross-platform checks (CRAN requires these)
```r
devtools::check_win_devel()      # emails you a result link
devtools::check_win_release()
rhub::rhub_check()               # Linux + macOS
```

### 4. Submit
```r
devtools::release()              # or upload tarball at https://cran.r-project.org/submit.html
```

---

## ⚠️ Known issues the check WILL flag (fix as they appear)

1. **Undefined globals** — many functions use base funcs without namespacing →
   "no visible global function/binding" NOTEs. Fix via `@importFrom` (step 1)
   and/or `utils::globalVariables()` for data-frame column names used in
   non-standard evaluation.
2. **Examples** — `\dontrun{}` is used heavily. CRAN dislikes hiding runnable
   code in `\dontrun{}`. For examples that are merely slow, switch to
   `\donttest{}`. Keep `\dontrun{}` ONLY for interactive
   (`drawExtent`/`drawPoly`/`readline`) or internet/file-writing examples.
   Provide at least one small, fast, runnable example per exported function.
3. **Internet access** (`weather.R`: NASA/WMO downloads) — must use `https`,
   wrap in `tryCatch()`, and fail gracefully (return `NULL`/message, not
   error) when offline. Never run in examples/tests without a guard.
4. **File writing** (`rsebWrite.R`, `weather.R` `write.csv`/`writeRaster`) —
   fine as functions (they write to a user-supplied `folder`), but their
   examples must write only to `tempdir()`.
5. **`T`/`F`** instead of `TRUE`/`FALSE` — tidy up (style NOTE in some checks).
6. **Spelling** — run `devtools::spell_check()`; add real terms to
   `inst/WORDLIST`.
7. **`\value`** — every exported function's `.Rd` must document a return value
   or you get a WARNING.

---

## 🔭 terra / sf migration (future-proofing — NOT a CRAN blocker today)

`raster` and `sp` are **still on CRAN and pass checks** — only `rgdal`,
`rgeos`, and `maptools` were retired (Oct 2023), and this package does **not**
use any of them. So you can submit on `raster` + `sp` now.

`sp` is in maintenance mode and `raster` is built on `terra`, so migrating is
worthwhile eventually. It is a large, pervasive rewrite (`raster()`→`rast()`,
`stack()`→`c()`/`rast()`, `getValues`→`values`, `crop`/`zonal`/`extract`,
`rasterToPoints`, `cellFromXY`, `CRS`/`proj4string`→`sf::st_crs`, `drawExtent`,
…) and should be done **with the test data in front of you** so numeric output
can be compared before/after. Recommended: submit now, migrate in 1.1-0.

## Submission summary

Resubmission of sebkc 1.0-6.

Addresses the review by Konstanze Lauseker: coldTs() and hotTs() no longer
set a fixed seed internally. Each function now takes a `seed` argument
(default NULL) that is passed to set.seed() only when the user supplies it,
so the random-number stream is left untouched by default and reproducibility
is under the user's control. The only remaining NOTE is the expected
"New submission" (the flagged words are model acronyms, author names and
"et al.").

## History

* 1.0-5: removed the bundled non-standard file 'CRAN-SUBMISSION' via
  .Rbuildignore.
* 1.0-4: removed an FAO URL that returned HTTP 404; cite FAO papers instead.
* 1.0-3: addressed the review by Konstanze Lauseker -- explained acronyms;
  T/F -> TRUE/FALSE; print() -> message(); removed options(warn=-1);
  restore par() and RNG state via on.exit(); fixed unexecutable example
  code and a length-1 bug in biomass(). Pure-numeric FAO-56 examples (ETo,
  ETohr) run under \donttest{}; heavy satellite / internet examples remain
  in \dontrun{}.

## Test environments

* local: Windows 10, R 4.6.0 -- 0 errors | 0 warnings | 0 notes

## R CMD check results

0 errors | 0 warnings | 1 note (New submission; installed size ~5.5 Mb due
to the bundled Landsat scene used by examples).

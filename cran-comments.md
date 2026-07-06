## Submission summary

Resubmission of sebkc 1.0-4.

This fixes the invalid URL reported by Uwe Ligges for 1.0-3:
the FAO crop-information page returned HTTP 404. The two URLs in kc.Rd /
cal.kc.Rd have been removed and replaced with a plain-text citation
(FAO Irrigation and Drainage Papers 33 and 56), so there is no URL to check.

## Earlier review (Konstanze Lauseker), addressed in 1.0-3 and retained

* Explained all acronyms in the Description.
* Replaced T/F with TRUE/FALSE; print() with message(); removed
  options(warn=-1); restore par() and the user's RNG state via on.exit();
  fixed unexecutable example code and a length-1 bug in biomass().
* Pure-numeric FAO-56 examples (ETo, ETohr) run under \donttest{}. The
  remaining examples stay in \dontrun{} because they need internet
  (weather) or process a full Landsat scene taking minutes each.

## Test environments

* local: Windows 10, R 4.6.0 -- 0 errors | 0 warnings | 0 notes

## R CMD check results

0 errors | 0 warnings | 1 note (New submission; installed size ~5.5 Mb due
to the bundled Landsat scene used by examples).

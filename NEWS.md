# tsissm 1.0.3

* Work around RTMB AD eigen failures on M1 macOS. CRAN M1MAC/R-devel checks fail 
while taping RTMB's AD overload for non-symmetric eigen decompositions, with an 
invalid `advector` error inside the complex eigenvector normalization path. 
This appears to be platform and toolchain sensitive: Apple silicon, especially M1 (but not M3), 
does not provide the same extended precision behavior as some x86 paths, and the 
reduced floating-point headroom can expose fragile complex eigen/AD calculations.
Added a lazy runtime capability check for RTMB AD eigen support. Normal platforms
continue to use RTMB's AD eigen constraints. If the capability check fails, the
package falls back only for eigenvalue constraints: constraint values are
computed with numeric `base::eigen()` and constraint Jacobians are computed by
central finite differences. The likelihood, scores, and other TMB/RTMB AD paths
remain unchanged.
* Exposed `tsissm_capabilities()` to report the detected RTMB AD eigen capability,
source of the decision, RTMB/TMB versions, and platform details.
* Documented the `tsissm.rtmb_ad_eigen` option under `estimate()` so users can force RTMB AD
eigen constraints or force the finite-difference fallback for diagnostics.
* Also avoided a fragile vignette conversion through `as_flextable()` by printing
the model summary directly.


# tsissm 1.0.2

* The SLSQP solver (from nloptr) sometimes will return NaN for the parameters during
the optimization phase. Added a stop clause else it may crash the process (due 
to RTMB eigen crashing).
* Added a solver control function (issm_control)

# tsissm 1.0.1

* Initial CRAN submission based on complete re-write of model on github which was
originally created on 2020-09-16.


# tsissm <img src="man/figures/logo.png" align="right" height="139" alt="" />

[![R-CMD-check](https://github.com/tsmodels/tsissm/actions/workflows/rcmdcheck.yaml/badge.svg)](https://github.com/tsmodels/tsissm/actions/workflows/rcmdcheck.yaml)
[![Last-changedate](https://img.shields.io/badge/last%20change-2025--07--12-yellowgreen.svg)](/commits/master)
[![packageversion](https://img.shields.io/badge/Package%20version-1.0.2-orange.svg?style=flat-square)](commits/master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/tsissm)](https://cran.r-project.org/package=tsissm)

# tsissm

Unobserved components model using the linear innovations state space
representation (single source of error).

Key features:

- Estimation using autodiff (TMB)
- Model selection and ensembling
- Regressors in the observation equation (non-time varying)
- Choice of distributions (Gaussian, Student and Johnson’s SU)
- Choice of dynamics in variance (constant or GARCH)

Methods for specification (`modelspec`), estimation (`estimate`),
summary with choice of vcov (`summary` and `vcov`), diagnostics
(`tsdiagnose`), online filtering (`tsfilter`), prediction (`predict`),
simulation (`simulate`), backtesting (`tsbacktest`) and profiling
(`tsprofile`). Additionally, model selection (`auto_select`) and
ensembling (`tsensemble`) is also included.

# tsissm 1.0.2

* The SLSQP solver (from nloptr) sometimes will return NaN for the parameters during
the optimization phase. Added a stop clause else it may crash the process (due 
to RTMB eigen crashing).
* Added a solver control function (issm_control)

# tsissm 1.0.1

* Initial CRAN submission based on complete re-write of model on github which was
originally created on 2020-09-16.

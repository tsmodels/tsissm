data("us_retail_sales")
suppressMessages(suppressWarnings(library(xts)))
suppressMessages(suppressWarnings(library(tsaux)))
y <- as.xts(us_retail_sales)

spec_constant_benchmark <- issm_modelspec(y, slope = TRUE, seasonal = TRUE, seasonal_frequency = 12, seasonal_harmonics = 5)
mod_constant_benchmark <- estimate(spec_constant_benchmark)
constant_benchmark <-  c("alpha" = 0.3761183584, "beta" = 0.0028909622, "gamma12.1" = -0.0005341079, 
                         "gamma12.2" = 0.0005066655, "sigma" = 12448.09, "LogLik" = -4317.608)

spec_dynamic_benchmark <- issm_modelspec(y, slope = TRUE, seasonal = TRUE, variance = "dynamic", seasonal_frequency = 12, seasonal_harmonics = 5)
mod_dynamic_benchmark <- estimate(spec_dynamic_benchmark)
dynamic_benchmark <-  c("alpha" = 0.2794620879, "beta" = 0.0009843874, "gamma12.1" = 0.0002238915, 
                        "gamma12.2" = -0.0009547473, "eta" = 0.1366928283, "delta" = 0.8237582308,
                        "LogLik" = -4270.23)

spec_bc_benchmark <- issm_modelspec(y, slope = TRUE, seasonal = TRUE, seasonal_frequency = 12, seasonal_harmonics = 5, lambda = 0.25)
mod_bc_benchmark <- estimate(spec_bc_benchmark)
bc_benchmark <-  c("alpha" = 0.400342866, "beta" = 0.003315568, "gamma12.1" = -0.006163750, 
                   "gamma12.2" = 0.001806082, "sigma" = 0.6742348, "LogLik" = -4186.793)

spec_auto <- issm_modelspec(y, auto = TRUE, slope = c(TRUE, FALSE), seasonal = TRUE, seasonal_frequency = 12, seasonal_harmonics = list(c(5)), top_n = 2)
mod_auto <- estimate(spec_auto)


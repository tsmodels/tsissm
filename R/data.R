#' Mauna Loa CO2 Daily Measurements
#'
#' The dataset represents the daily average of the number of carbon dioxide 
#' molecules in a given number of molecules of air, after removal of water 
#' vapor (molefrac) from the Mauna Loa observatory, spanning the period
#' 1974-05-19 to 2025-02-25.
#' @format ## `maunaloa`
#' A data.frame containing 15571 observations
#' \describe{
#'   \item{date}{The YYYY-MM-DD date}
#'   \item{co2}{The CO2 measurements}
#' }
#' @source \url{https://gml.noaa.gov/ccgg/trends/data.html}
"maunaloa"


#' Advance Retail Sales: Retail Trade (US)
#'
#' The dataset represents the monthly, non-seasonally adjusted advance estimate
#' of sales from the retail trade based on data from a subsample of firms 
#' from the larger Monthly Retail Trade Survey.
#' @format ## `us_retail_sales`
#' A data.frame containing 398 observations
#' \describe{
#'   \item{date}{The YYYY-MM-DD date}
#'   \item{co2}{The sales value in millions of dollars.}
#' }
#' @source U.S. Census Bureau, Advance Retail Sales: Retail Trade [RSXFSN]. Retrieved from FRED, Federal Reserve Bank of St. Louis: \url{https://fred.stlouisfed.org/series/RSXFSN} on March 27, 2025.
"us_retail_sales"


#' SPY ETF Adjusted Close
#'
#' The adjusted daily closing price of the SPY ETF.
#' 
#' @format ## `spy`
#' A data.frame containing 7597 observations
#' \describe{
#'   \item{date}{The YYYY-MM-DD date}
#'   \item{spy}{The adjusted stock closing price.}
#' }
#' @source From Yahoo Finance.
"spy"
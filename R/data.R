#' US Treasury bill and bond rates (term-structure example)
#'
#' Monthly 3-month Treasury bill rate and 10-year Treasury constant-maturity
#' rate, the data of the applied example of the ATSE cointegration chapter
#' (the expectations hypothesis of the term structure). The sample starts
#' when both series exist (April 1953).
#'
#' @format A data frame with one row per month and 3 columns:
#' \describe{
#'   \item{date}{first day of the month (`Date`).}
#'   \item{tb3ms}{3-Month Treasury Bill Secondary Market Rate, Discount Basis,
#'     percent, not seasonally adjusted (FRED `TB3MS`).}
#'   \item{gs10}{Market Yield on U.S. Treasury Securities at 10-Year Constant
#'     Maturity, Quoted on an Investment Basis, percent, not seasonally
#'     adjusted (FRED `GS10`).}
#' }
#' The download date is in `attr(term_structure, "retrieved")`.
#' @source Board of Governors of the Federal Reserve System (US), H.15
#'   Selected Interest Rates, retrieved from FRED, Federal Reserve Bank of
#'   St. Louis, <https://fred.stlouisfed.org/>. Script: `data-raw/term_structure.R`.
#' @examples
#' head(term_structure)
#' eg_test(gs10 ~ tb3ms, term_structure)
"term_structure"

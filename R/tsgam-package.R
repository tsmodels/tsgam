#' @keywords internal
#' @import methods
#' @import tsmethods
#' @import data.table
#' @import mgcv
#' @importFrom stats as.formula gaussian quantile na.omit predict residuals coef fitted AIC sd sigma cor
#' @importFrom gratia predicted_samples
#' @importFrom tsaux smape mape bias crps msis mis mslre sampling_frequency
#' @importFrom future.apply future_lapply
#' @importFrom rlang f_rhs f_lhs
#' @importFrom copula normalCopula tCopula P2p rCopula fitCopula pobs
#' @importFrom future %<-%
#' @importFrom tsdistributions distribution_modelspec rdist qdist
#' @importFrom progressr handlers progressor
#' @importFrom zoo index as.zoo zoo coredata
#' @importFrom fasttime fastDate fastPOSIXct
#' @importFrom xts xts as.xts is.xts tzone
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

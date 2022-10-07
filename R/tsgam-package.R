#' @keywords internal
#' @import methods
#' @import tsmethods
#' @import data.table
#' @import mgcv
#' @importFrom stats as.formula gaussian quantile na.omit predict residuals
#' @importFrom gratia predicted_samples
#' @importFrom tsaux smape mape bias crps msis mis sampling_frequency
#' @importFrom future.apply future_lapply
#' @importFrom future %<-%
#' @importFrom tsdistributions distribution_modelspec rdist
#' @importFrom progressr handlers progressor
#' @importFrom zoo index as.zoo zoo coredata
#' @importFrom xts xts as.xts is.xts tzone
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

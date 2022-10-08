#' Model Estimation Summary
#'
#' @description Summary method for class \dQuote{gam.estimate}
#' @param object an object of class \dQuote{gam.estimate}.
#' @param ... additional arguments passed through to the mgcv gam summary method.
#' @return A list of summary information (see mgcv summary method).
#' @aliases summary
#' @method summary gam.estimate
#' @rdname summary
#' @export
#'
#'
summary.gam.estimate <- function(object, ...)
{
  return(summary(object$model, ...))
}

#' Model Fitted Values
#'
#' @description Extract the fitted values from an estimated model.
#' @param object an object of class \dQuote{gam.estimate}.
#' @param ... not currently used.
#' @return An xts vector
#' @aliases fitted
#' @method fitted gam.estimate
#' @rdname fitted
#' @export
#'
#'
fitted.gam.estimate <- function(object, ...)
{
  fit_values <- object$model$y - as.numeric(residuals(object$model, type = "response"))
  out <- xts(fit_values, object$spec$clean_time_index)
  colnames(out) <- "fitted"
  return(out)
}

#' Model Residuals
#'
#' @description Extract the residual values from an estimated model.
#' @param object an object of class \dQuote{gam.estimate}.
#' @param type the of residuals with \dQuote{response} being the default. Other
#' options see the mgcv documentation.
#' @param ... not currently used.
#' @return An xts vector
#' @aliases residuals
#' @method residuals gam.estimate
#' @rdname residuals
#' @export
#'
#'
residuals.gam.estimate <- function(object, type = "response", ...)
{
  resid_values <- as.numeric(residuals(object$model, type = type))
  out <- xts(resid_values, object$spec$clean_time_index)
  colnames(out) <- "residuals"
  return(out)
}

#' Performance Metrics
#'
#' @description Performance metrics from an estimated or predicted gam model.
#' @param object an object of class \dQuote{gam.estimate} or \dQuote{gam.predict}
#' @param actual the actual data matched to the dates of the forecasts.
#' @param alpha the coverage level for distributional forecast metrics.
#' @param ... not currently used.
#' @aliases tsmetrics
#' @method tsmetrics gam.predict
#' @rdname tsmetrics
#' @export
#'
#'
tsmetrics.gam.predict <- function(object, actual, alpha = 0.1, ...)
{
  n <- NCOL(object$distribution)
  if (NROW(actual) != n) stop("\nactual length not equal to forecast output length")
  m_mape <- mape(actual, object$mean)
  m_smape <- smape(actual, object$mean)
  m_bias <- bias(actual, object$mean)
  m_mslre <- mslre(actual, object$mean)
  if (!is.null(alpha)) {
    m_mis <- mis(actual, lower = apply(object$distribution, 2, quantile, alpha/2), upper = apply(object$distribution, 2, quantile, 1 - alpha/2), alpha = alpha)
  } else {
    m_mis <- as.numeric(NA)
  }
  m_crps <- crps(actual, object$distribution)
  return(data.frame("MAPE" = m_mape, "SMAPE" = m_smape,
             "MSLRE" = m_mslre, "BIAS" = m_bias,
             "MIS" = m_mis, "CRPS" = m_crps))
}

#' @method tsmetrics gam.estimate
#' @rdname tsmetrics
#' @export
#'
#'
tsmetrics.gam.estimate = function(object, ...)
{
  # residuals diagnostics
  fx <- fitted(object)
  ac <- xts(object$model$y, object$spec$clean_time_index)
  acfx <- na.omit(cbind(ac, fx))
  actual <- as.numeric(acfx[,1])
  fted <- as.numeric(acfx[,2])
  m_mape <- mape(actual, fted)
  m_smape <- smape(actual, fted)
  m_mslre  <- mslre(actual, fted)
  m_bias <- bias(actual, fted)
  m_sum <- summary(object)
  m_deviance <- m_sum$dev.expl

  return(data.frame("n" = m_sum$n,
                    "no.pars" = m_sum$np,
                    "Deviance_Explained" = m_deviance,
                    "AIC" = AIC(object$model),
                    "MAPE" = m_mape,
                    "SMAPE" = m_smape,
                    "MSLRE" = m_mslre,
                    "BIAS" = m_bias))
}

decompose_model <- function(object, newdata = NULL, type = "estimate")
{
  if (type == "estimate") {
    terms <- predict(object$model, type = "lpmatrix")
  } else {
    terms <- predict(object$model, newdata = newdata, type = "lpmatrix")
  }
  coefficients <- coef(object$model)
  cnames <- colnames(terms)
  sumr <- summary(object$model)
  unique_terms <- rownames(sumr$s.table)
  if (any(cnames %in% "(Intercept)")) unique_terms <- c("(Intercept)", unique_terms)
  output_names <- unique_terms
  unique_terms <- gsub("\\(","\\\\(", unique_terms)
  unique_terms <- gsub("\\)","\\\\)", unique_terms)
  decomp <- do.call(cbind, lapply(1:length(unique_terms), function(i){
    use <- which(grepl(unique_terms[i], cnames, perl = TRUE))
    x <- terms[,use, drop = FALSE] %*% coefficients[use]
    return(x)
  }))
  colnames(decomp) <- output_names
  return(decomp)
}

#' Model Decomposition
#'
#' @description Decomposes the estimated model or prediction into its component
#' parts (aggregated basis terms).
#' @param object an object of class \dQuote{gam.estimate} or \dQuote{gam.predict}
#' @param ... not currently used.
#' @return An xts matrix of the decomposed additive components.
#' @aliases tsdecompose
#' @method tsdecompose gam.estimate
#' @rdname tsdecompose
#' @export
#'
#'
tsdecompose.gam.estimate <- function(object, ...)
{
  decomp <- decompose_model(object, type = "estimate")
  decomp <- xts(decomp, object$spec$clean_time_index)
  return(decomp)
}

#' @method tsdecompose gam.predict
#' @rdname tsdecompose
#' @export
#'
tsdecompose.gam.predict <- function(object,  ...)
{
  return(object$decomposition)
}

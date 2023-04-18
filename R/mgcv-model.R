#' GAM Model specification
#'
#' @description Sets up a gam model via formula.
#' @param formula a valid gam formula.
#' @param data an xts matrix of the data.
#' @param family a valid gam distribution famiily
#' @param ... additional arguments (none currently supported).
#' @return An object of class \dQuote{gam.spec}.
#' @note This is a wrapper to the gam function from the mgcv package.
#' @aliases gam_modelspec
#' @rdname gam_modelspec
#' @export
#'
gam_modelspec <- function(formula, data, family = gaussian(link = "identity"), ...)
{
  if (!is.xts(data)) {
    stop("data must be an xts object")
  }
  if (!is(formula,"formula")) {
    formula <- as.formula(formula)
  }
  # check formula against data
  formula_vars <- all.vars(formula)
  data_vars <- colnames(data)
  data <- data[,formula_vars]
  check <- all.equal(formula_vars, colnames(data))
  if (!check) stop("\ndata variables do not match formula variables/")
  time_zone <- tzone(index(data))
  spec <- list()
  spec$formula <- formula
  spec$data <- as.data.frame(data, row.names = FALSE)
  spec$time_index <- index(data)
  spec$clean_time_index <- index(na.omit(data))
  spec$time_class <- attr(sampling_frequency(index(data)),"date_class")
  spec$time_zone <- time_zone
  spec$family <- family
  class(spec) <- "gam.spec"
  return(spec)
}


#' Model Estimation
#'
#' @description Estimates a model given a specification object.
#' @param object an object of class \dQuote{gam.spec}.
#' @param ... additional arguments passed to the gam function of the mgcv package.
#' @details The object class retains the output model from mgcv in the model slot
#' of the returned object list so can be dispatched to other methods created for
#' mgcv gam objects.
#' @return An object of class \dQuote{gam.estimate}
#' @aliases estimate
#' @method estimate gam.spec
#' @rdname estimate
#' @export
#'
#'
estimate.gam.spec <- function(object, ...)
{
  mod <- gam(formula = object$formula, family = object$family, data = object$data, ...)
  extra <- list()
  extra$time_index <- object$time_index
  extra$time_class <- object$time_class
  extra$clean_time_index <- object$clean_time_index
  extra$time_zone <- object$time_zone
  out <- list(model = mod, spec = extra)
  class(out) <- "gam.estimate"
  return(out)
}

#' Model Prediction
#'
#' @description Prediction function for class \dQuote{gam.estimate}.
#' @param object an object of class \dQuote{gam.estimate}.
#' @param newdata an xts matrix of external regressors in the forecast horizon.
#' @param nsim The number of simulations to use for generating the simulated predictive distribution.
#' @param distribution A valid distribution from the tsdistributions package to be used to fit the
#' response residuals and then proxy for the predictive distribution. Only valid if the estimated
#' model used the gaussian family, else will throw an error.
#' @param innov an optional matrix of innovations (see innov_type for admissible types),
#' of dimensions nsim x horizon.
#' @param innov_type if \\sQuote{innov} is not NULL, then this denotes the type of values
#' passed, with \\dQuote{q} denoting quantile probabilities (default and
#' backwards compatible), \\dQuote{z} for standardized errors and \\dQuote{r} for
#' raw innovations (non standardized). See details.
#' @param tabular whether to return the data in long format as a data.table instead
#' of an object of tsmodels.predict
#' @param ... not currently used.
#' @return An object which inherits class \dQuote{tsmodel.predict} with slots for
#' the simulated predictive distribution, the original series (as a
#' zoo object), the original specification object and the mean forecast. If tabular
#' is TRUE, will return a data.table object of the predictions in long format with
#' some additional columns.
#' Using a gaussian distribution for estimation followed by a custom distribution for fitting
#' the model residuals in the prediction step to proxy for the uncertainty is a 2 step option
#' which may yield superior results to using something like family scat in the mgcv model as this
#' has proven to be problematic in certain cases.
#' @details
#' The \dQuote{distribution} and \dQuote{innovations} arguments allow for significant customization
#' of the simulated distribution. When \dQuote{innov_type} is type q, then the values
#' are transformed back to normally distributed predictions using the normal quantile function
#' with mean the prediction vector and standard deviation calculated from the model residuals.
#' If distribution is also not NULL, then instead of using the normal quantile function, a model
#' if first fitted to the residuals given the distribution chosen and the quantile function of
#' that distribution is used to transform the uniform value.
#' When \dQuote{innov_type} is either z or r, then the values are simply rescaled using the
#' standard deviation of the residuals (for type z), and re-centered by the prediction means
#' (for type z and r).
#' @aliases predict
#' @method predict gam.estimate
#' @rdname predict
#' @export
#'
#'
predict.gam.estimate <- function(object, newdata, nsim = 9000, tabular = FALSE, distribution = NULL, innov = NULL, innov_type = "q", ...)
{
  # setup data.table variable names to avoid notes in checking
  prediction <- parameter <- forecast_date <- estimation_date <- draw <- value <- NULL
  `.` <- list
  newx <- as.data.frame(newdata)
  if (any(is.na(newdata))) stop("\nno NAs allowed in data")
  mean_prediction <- predict(object$model, newdata = newx)
  if (is.matrix(mean_prediction)) mean_prediction <- prediction[,1]
  # create the decomposition
  decomp <- decompose_model(object, newdata = newx, type = "predict")
  decomp <- xts(decomp, index(newdata))
  if (!is.null(innov) & innov_type != "q") distribution <- NULL
  if (!is.null(distribution)) {
    if (object$model$family$family != "gaussian") stop("\ndistribution option only available for models estimated using gaussian family")
    # estimate distribution on residuals of model
    # predict given mu and distributional parameters
    spec <- distribution_modelspec(residuals(object$model, type = "response"), distribution = distribution)
    distribution_fit <- estimate(spec)
    p_matrix <- distribution_fit$spec$parmatrix
    sigma <- p_matrix[parameter == "sigma"]$value
    skew <- p_matrix[parameter == "skew"]$value
    shape <- p_matrix[parameter == "shape"]$value
    lambda <- p_matrix[parameter == "lambda"]$value
    simulated_draws <- do.call(cbind, lapply(1:length(mean_prediction), function(i) {
      rdist(distribution = distribution, n = nsim, mu = mean_prediction[i], sigma = sigma,
                             skew = skew, shape = shape, lambda = lambda)
    }))
    if (tabular) {
      simulated_draws <- as.data.table(simulated_draws)
      simulated_draws[,draw := 1:.N]
      simulated_draws <- melt(simulated_draws, id.vars = "draw", measure.vars = 1:(ncol(simulated_draws) - 1), variable.name = "forecast_date", value.name = "value")
      if (object$spec$time_class == "POSIXct") {
        fun <- as.POSIXct
      } else {
        fun <- as.Date
      }
      simulated_draws[,forecast_date := fun(forecast_date, tz = object$spec$time_zone)]
      e_date <- max(object$spec$clean_time_index)
      simulated_draws[,estimation_date := e_date]
      simulated_draws <- simulated_draws[,.(estimation_date, draw, forecast_date, value)]
    } else {
      simulated_draws <- as.matrix(simulated_draws)
      colnames(simulated_draws) <- as.character(index(newdata))
      class(simulated_draws) <- "tsmodel.distribution"
      attr(simulated_draws, "date_class") <- object$time_class
      simulated_draws <- list(original_series = zoo(object$model$y, object$spec$clean_time_index),
                              distribution = simulated_draws,
                              mean = zoo(mean_prediction, index(newdata)), decomposition = decomp)
      class(simulated_draws) <- c("gam.predict","tsmodel.predict")
      return(simulated_draws)
    }
  } else {
    if (!is.null(innov)) {
      h <- NROW(newdata)
      innov_type <- match.arg(innov_type, c("q", "z", "r"))
      innov <- as.matrix(innov)
      if (NROW(innov) != nsim) {
        stop("\nnrow of innov must be nsim")
      }
      if (NCOL(innov) != h) {
        stop("\nncol of innov must be NROW(newdata)")
      }
      # check that the innovations are uniform samples (from a copula)
      if (innov_type == "q") {
        if (any(innov < 0 | innov > 1)) {
          stop("\ninnov must be >0 and <1 (uniform samples) for innov_type = 'q'")
        }
        if (any(innov == 0))
          innov[which(innov == 0)] <- 1e-12
        if (any(innov == 1))
          innov[which(innov == 1)] <- (1 - 1e-12)
      }
      innov <- matrix(innov, nsim, h)
      sig <- sd(object$model$residuals)
      if (innov_type == "q") {
        if (!is.null(distribution)) {
          spec <- distribution_modelspec(residuals(object$model, type = "response"), distribution = distribution)
          distribution_fit <- estimate(spec)
          p_matrix <- distribution_fit$spec$parmatrix
          sigma <- p_matrix[parameter == "sigma"]$value
          skew <- p_matrix[parameter == "skew"]$value
          shape <- p_matrix[parameter == "shape"]$value
          lambda <- p_matrix[parameter == "lambda"]$value
          simulated_draws <- do.call(cbind, lapply(1:length(mean_prediction), function(i) {
            qdist(distribution = distribution, p = innov[,i], mu = mean_prediction[i], sigma = sigma,
                  skew = skew, shape = shape, lambda = lambda)
          }))
        } else {
          simulated_draws <- do.call(cbind, lapply(1:length(mean_prediction), function(i) {
            qdist(distribution = "norm", p = innov[,i], mu = mean_prediction[i], sigma = sig)
            }))
        }
      } else {
        if (innov_type == "z") {
          simulated_draws <- sweep(innov, 2, sig, "*")
          simulated_draws <- sweep(simulated_draws, 2, mean_prediction, "+")
        } else {
          simulated_draws <- sweep(innov, 2, mean_prediction, "+")
        }
      }
      if (tabular) {
        simulated_draws <- as.data.table(simulated_draws)
        simulated_draws[,draw := 1:.N]
        simulated_draws <- melt(simulated_draws, id.vars = "draw", measure.vars = 1:(ncol(simulated_draws) - 1), variable.name = "forecast_date", value.name = "value")
        if (object$spec$time_class == "POSIXct") {
          fun <- as.POSIXct
        } else {
          fun <- as.Date
        }
        simulated_draws[,forecast_date := fun(forecast_date, tz = object$spec$time_zone)]
        e_date <- max(object$spec$clean_time_index)
        simulated_draws[,estimation_date := e_date]
        simulated_draws <- simulated_draws[,.(estimation_date, draw, forecast_date, value)]
      } else {
        simulated_draws <- as.matrix(simulated_draws)
        colnames(simulated_draws) <- as.character(index(newdata))
        class(simulated_draws) <- "tsmodel.distribution"
        attr(simulated_draws, "date_class") <- object$time_class
        simulated_draws <- list(original_series = zoo(object$model$y, object$spec$clean_time_index),
                                distribution = simulated_draws,
                                mean = zoo(mean_prediction, index(newdata)), decomposition = decomp)
        class(simulated_draws) <- c("gam.predict","tsmodel.predict")
        return(simulated_draws)
      }
    } else {
      simulated_draws <- as.data.table(predicted_samples(object$model, n = nsim, newdata = newx))
      simulated_draws <- dcast(simulated_draws, draw~row, value.var = "response")
      simulated_draws <- simulated_draws[order(draw)]
      if (tabular) {
        setcolorder(simulated_draws, c("draw",paste0(1:(ncol(simulated_draws) - 1))))
        colnames(simulated_draws) <- c("draw", as.character(index(newdata)))
        simulated_draws <- melt(simulated_draws, id.vars = "draw", measure.vars = 2:ncol(simulated_draws), variable.name = "forecast_date", value.name = "value")
        if (object$spec$time_class == "POSIXct") {
          fun <- as.POSIXct
        } else {
          fun <- as.Date
        }
        simulated_draws[,forecast_date := fun(forecast_date, tz = object$spec$time_zone)]
        e_date <- max(object$spec$clean_time_index)
        simulated_draws[,estimation_date := e_date]
        simulated_draws <- simulated_draws[,.(estimation_date, draw, forecast_date, value)]
        return(simulated_draws)
      } else {
        simulated_draws <- simulated_draws[,draw := NULL]
        setcolorder(simulated_draws, paste0(1:ncol(simulated_draws)))
        simulated_draws <- as.matrix(simulated_draws)
        colnames(simulated_draws) <- as.character(index(newdata))
        class(simulated_draws) <- "tsmodel.distribution"
        attr(simulated_draws, "date_class") <- object$time_class
        simulated_draws <- list(original_series = zoo(object$model$y, object$spec$clean_time_index),
                                distribution = simulated_draws, mean = zoo(mean_prediction, index(newdata)),decomposition = decomp)
        class(simulated_draws) <- c("gam.predict","tsmodel.predict")
        return(simulated_draws)
      }
    }
  }
}


#' GAM Model Training Data Setup
#'
#' @description Sets up a gam data object for backtesting.
#' @param formula a valid gam formula.
#' @param family a valid gam distribution family
#' @param data an xts data matrix
#' @param estimation_dates a vector of estimation dates
#' @param prediction_dates a list of length(estimation_dates) with the
#' dates to use for forecasting after each estimation date.
#' @param validate whether to perform validity tests on the prediction_dates
#' so that they are never less than or equal to estimation_dates, and that both
#' estimation_dates and prediction_dates actually exist in the supplied data. The
#' first one will generate a wwrning, and the last 2 an error.
#' @return an object of class \dQuote{gam.trainspec}
#' @aliases gam_trainspec
#' @rdname gam_trainspec
#' @export
#'
#'
gam_trainspec <- function(formula, family = gaussian(link = "identity"), data, estimation_dates, prediction_dates, validate = TRUE)
{
  if (!is.xts(data)) {
    stop("data must be an xts object")
  }
  if (!is(formula,"formula")) {
    formula <- as.formula(formula)
  }
  if (!is.list(prediction_dates)) stop("\nprediction_dates must be a list")
  n <- length(estimation_dates)
  if (length(prediction_dates) != n) stop("\nlength of prediction_dates list must be equal to length(estimation_dates)")
  if (validate) {
    check <- sapply(1:n, function(i){
      all(prediction_dates[[i]] > estimation_dates[i])
    })
    if (any(!check)) warning("\nvalidation failed for prediction_dates>estimation_dates.")
    data_dates <- index(data)
    check <- sapply(1:n, function(i){
      all(estimation_dates[i] %in% data_dates)
    })
    if (any(!check)) stop("\nvalidation failed for estimation_dates in data dates.")
    check <- sapply(1:n, function(i){
      all(prediction_dates[[i]] %in% data_dates)
    })
    if (any(!check)) stop("\nvalidation failed for prediction_dates in data dates.")
  }
  # check formula against data
  formula_vars <- all.vars(formula)
  data_vars <- colnames(data)
  data <- data[,formula_vars]
  check <- all.equal(formula_vars, colnames(data))
  if (!check) stop("\ndata variables do not match formula variables/")
  spec <- list()
  spec$formula <- formula
  spec$family <- family
  spec$data <- data
  spec$estimation_dates <- estimation_dates
  spec$prediction_dates <- prediction_dates
  class(spec) <- "gam.trainspec"
  return(spec)
}

#' Walk Forward Model Backtest
#'
#' @description Generates an expanding window walk forward backtest.
#' @param object an object of class \dQuote{gam.train}.
#' @param alpha optional numeric vector of coverage rates for which to calculate
#' the quantiles.
#' @param distribution A valid distribution from the tsdistributions package to be used to fit the
#' response residuals and then proxy for the predictive distribution. Only valid if the estimated
#' model used the gaussian family, else will throw an error.
#' @param tabular whether to return the full simulated distributiion in long format
#' as a data.table instead of summarized distribution with quantiles.
#' @param nsim The number of simulations to use for generating the simulated predictive distribution.
#' @param trace whether to show the progress bar. The user is expected to have
#' set up appropriate handlers for this using the \dQuote{progressr} package.
#' @param ... additional arguments passed to the \dQuote{gam} function during estimation.
#' @return A list with the following data.tables:
#' \itemize{
#' \item prediction : the backtest table with forecasts and actuals
#' \item metrics: a summary performance table showing metrics by
#' forecast horizon (MAPE, MSLRE, BIAS and MIS if alpha was not NULL).
#' }
#' @note The function can use parallel functionality as long as the user has
#' set up a \code{\link[future]{plan}} using the future package.
#' @aliases tsbacktest
#' @method tsbacktest gam.trainspec
#' @rdname tsbacktest
#' @export
#'
#'
tsbacktest.gam.trainspec <- function(object, alpha = NULL, trace = FALSE, nsim = 5000, distribution = NULL, tabular = FALSE, ...)
{
  if (!is.null(alpha)) {
    if (any(alpha <= 0)) {
      stop("\nalpha must be strictly positive")
    }
    if (any(alpha >= 1)) {
      stop("\nalpha must be less than 1")
    }
    quantiles <- as.vector(sapply(1:length(alpha), function(k) c(alpha[k]/2, 1 - alpha[k]/2)))
  } else {
    quantiles <- NULL
  }
  estimation_dates <- object$estimation_dates
  prediction_dates <- object$prediction_dates
  formula <- object$formula
  data <- object$data
  n <- length(estimation_dates)
  if (trace) {
    prog_trace <- progressor(n)
  }
  extra_args <- list(...)
  b %<-% future_lapply(1:n, function(i) {
    if (trace) prog_trace()
    data_train <- data[paste0("/", estimation_dates[i])]
    data_predict <- data[prediction_dates[[i]]]
    spec <- gam_modelspec(object$formula, data_train, family = object$family)
    mod <- do.call(estimate, args = list(object = spec, `...` = extra_args), quote = T)
    if (tabular) {
      out <- predict(mod, newdata = data_predict, distribution = distribution, nsim = nsim, tabular = TRUE)
    } else {
      p <- predict(mod, newdata = data_predict, distribution = distribution, nsim = nsim, tabular = FALSE)
      if (!is.null(quantiles)) {
        qp <- apply(p$distribution, 2, quantile, quantiles)
        if (length(quantiles) == 1) {
          qp <- matrix(qp, ncol = 1)
        } else{
          qp <- t(qp)
        }
        colnames(qp) <- paste0("P", round(quantiles*100, 2))
      }
      target <- all.vars(mod$model$formula)[1]
      out <- data.table("estimation_date" = estimation_dates[i],
                        "forecast_dates" = index(data_predict),
                        "forecast" = as.numeric(p$mean), "actual" = as.numeric(data_predict[,target]))
      if (!is.null(quantiles)) out <- cbind(out, qp)
    }
    # add crps and msis
    return(out)
  }, future.packages = c("tsmethods","tsaux","mgcv","xts","gratia","data.table"), future.seed = TRUE)
  b <- rbindlist(b)
  if (tabular) {
    return(b)
  } else {
    data_name <- "y"
    actual <- NULL
    forecast <- NULL
    metrics <- b[,list(variable = data_name, SMAPE = smape(actual, forecast), BIAS = bias(actual, forecast),
                       n = .N), by = "estimation_date"]
    if (!is.null(alpha)) {
      q_names <- matrix(paste0("P", round(quantiles*100, 2)), ncol = 2, byrow = TRUE)
      q <- do.call(cbind, lapply(1:length(alpha), function(i){
        b[,list(mis = mis(actual, get(q_names[i,1]), get(q_names[i,2]), alpha[i])), by = "estimation_date"]
      }))
      q <- q[,which(grepl("mis",colnames(q))), with = FALSE]
      colnames(q) <- paste0("MIS[",alpha,"]")
      metrics <- cbind(metrics, q)
    }
    return(list(prediction = b, metrics = metrics))
  }
}

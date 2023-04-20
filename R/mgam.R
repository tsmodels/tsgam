#' Multivariate GAM Model specification
#'
#' @description Sets up an mgam model via formula.
#' @param formula a valid gam formula.
#' @param data a data.table object in long format of the data, with no missing values.
#' @param series_id the column name which corresponds to the unique identifier
#' of the series.
#' @param index_id the column name which corresponds to the unique date/time
#' index of the series.
#' @param family a valid gam distribution family.
#' @param pooled whether to estimate a pooled model.
#' @param ... additional arguments (none currently supported).
#' @return An object of class \dQuote{gam.spec}.
#' @details
#' The multivariate GAM model described here is at a minimum a wrapper for multiple
#' estimation and prediction of a set of series in a convenient wrapper, with the
#' caveat that they share a common formula. However, it also allows pooled estimation
#' as well as forecasts with coupled uncertainty via the use of a copula.
#' @aliases mgam_modelspec
#' @rdname mgam_modelspec
#' @export
#'
mgam_modelspec <- function(formula, data, series_id = NULL, index_id = NULL, family = gaussian(link = "identity"), pooled = FALSE, ...)
{
  if (!inherits(data,"data.table")) stop("\ndata must be a data.table object")
  if (is.null(series_id)) stop("\nseries_id cannot be NULL")
  if (is.null(index_id)) stop("\nindex_id cannot be NULL")
  if (any(is.na(data))) {
    warning("\nmissing data found. Removing...")
    data <- na.omit(data)
  }
  if (family$family == "gaulss") {
    # formula must be list of length 2
    if (!is(formula, "list")) stop("\nfor the gaulss family, formula must be a list of length 2.")
    if (!is(formula[[1]],"formula")) formula[[1]] <- as.formula(formula[[1]])
    if (!is(formula[[2]],"formula")) formula[[2]] <- as.formula(formula[[2]])
    formula_vars <- unique(as.vector(unlist(sapply(formula, function(x) all.vars(f_rhs(x))))))
    response_var <- all.vars(f_lhs(formula[[1]]))
  } else if (family$family == "gaussian") {
    if (!is(formula,"formula")) {
      formula <- as.formula(formula)
    }
    formula_vars <- all.vars(f_rhs(formula))
    response_var <- all.vars(f_lhs(formula))
  }
  col_names <- colnames(data)
  series <- unique(data[[series_id]])
  # class of date column
  if (length(series) == 1) stop("\nnumber of unique ids but be > 1.")
  # check formula against data
  data_vars <- colnames(data)
  check <- all(c(response_var,formula_vars) %in% data_vars)
  if (!check) stop("\ndata variables do not match formula variables/")
  # remove any NA's
  data <- na.omit(data[,c(series_id, index_id, response_var, formula_vars), with = FALSE])
  data_index <- unique(data[[index_id]])
  time_zone <- tzone(data_index)
  spec <- list()
  spec$formula <- formula
  spec$pooled <- pooled
  spec$data <- data
  spec$series <- series
  spec$series_id <- series_id
  spec$index_id <- index_id
  spec$time_class <- attr(sampling_frequency(data_index),"date_class")
  spec$time_zone <- time_zone
  spec$response_var <- response_var
  spec$formula_vars <- formula_vars
  spec$family <- family
  class(spec) <- "mgam.spec"
  return(spec)
}

#' Model Estimation
#'
#' @description Estimates a model given a specification object.
#' @param object an object of class \dQuote{mgam.spec}.
#' @param use_bam whether to use \code{\link[mgcv]{bam}} for very large datasets
#' (in the pooled case only).
#' @param ... additional arguments passed to the gam function of the mgcv package, or
#' the bam function if \dQuote{use_bam} is true.
#' @details If the model is pooled, the object class retains the output model from
#' mgcv in the model slot of the returned object list so can be dispatched to
#' other methods created for mgcv gam objects. In the case of a non-pooled model,
#' the output id a list of models, one for each series.
#' For the non pooled model, the function makes use of the future package for
#' parallel estimation across the series. However, care should be taken in this
#' case since mgcv already makes use of openMP so the user may need to control
#' the number of threads already committed.
#' @return An object of class \dQuote{mgam.estimate}
#' @aliases estimate
#' @method estimate mgam.spec
#' @rdname estimate
#' @export
#'
#'
estimate.mgam.spec <- function(object, use_bam = FALSE, ...)
{
  .fitted <- .sigma <- .residuals <- .mu <- .sigma <- NULL
  pooled <- object$pooled
  if (pooled) {
    .x <- object$data
    if (use_bam) {
      mod <- bam(formula = object$formula, family = object$family, data = .x, ...)
    } else {
      mod <- gam(formula = object$formula, family = object$family, data = .x, ...)
    }
    .fit <- fitted(mod, type = "response")
    if (is.matrix(.fit)) {
      .mu <- .fit[,1]
      .sig <- 1/.fit[,2]
    } else {
      .mu <- .fit
      .sig <- sd(residuals(mod, type = "response"))
    }
    .x[,.fitted := .mu]
    .x[,.sigma := .sig]
    .x[,.residuals := as.numeric(residuals(mod, type = "response"))]
    spec <- object
    spec$data <- .x
    out <- list(model = mod, spec = spec)
    class(out) <- "mgam.estimate"
    return(out)
  } else {
    n <- length(object$series)
    .s <- split(object$data, by = c(object$series_id))
    mod %<-% future_lapply(1:n, function(i){
      .mod <- gam(formula = object$formula, family = object$family, data = .s[[i]], ...)
      return(.mod)
    }, future.packages = c("tsmethods","tsgam","mgcv","data.table"), future.seed = TRUE)
    mod <- eval(mod)
    .mu <- as.vector(sapply(mod, function(x) {
      .fit <- fitted(x, type = "response")
      if (is.matrix(.fit)) {
        .fit[,1]
      } else {
        .fit
      }
    }))
    .sig <- as.vector(sapply(mod, function(x) {
      .fit <- fitted(x, type = "response")
      if (is.matrix(.fit)) {
        1/.fit[,2]
      } else {
        rep(sd(residuals(x, type = "response")), length(.fit))
      }
    }))
    .res <- as.vector(sapply(mod, function(x) as.numeric(residuals(x, type = "response"))))
    .x <- object$data
    .x[,.fitted := .mu]
    .x[,.sigma := .sig]
    .x[,.residuals := .res]
    spec <- object
    spec$data <- .x
    spec$pooled <- FALSE
    out <- list(model = mod, spec = spec)
    class(out) <- "mgam.estimate"
    return(out)
  }
}

#' Extract Model Fitted Values
#'
#' @description Extract the fitted values of the estimated model.
#' @param object an object of class \dQuote{mgam.estimate}.
#' @param ... not currently used.
#' @return A data.table with the index_id, series_id and response fitted values.
#' @aliases fitted
#' @method fitted mgam.estimate
#' @rdname fitted
#' @export
#'
#'
fitted.mgam.estimate <- function(object, ...)
{
  out <- object$spec$data[,c(object$spec$series_id, object$spec$index_id, ".fitted"), with = FALSE]
  return(out)
}

#' Extract Model Residuals
#'
#' @description Extract the residuals of the estimated model.
#' @param object an object of class \dQuote{mgam.estimate}.
#' @param ... not currently used.
#' @return A data.table with the index_id, series_id and response residuals values.
#' @aliases residuals
#' @method residuals mgam.estimate
#' @rdname residuals
#' @export
#'
#'
residuals.mgam.estimate <- function(object, ...)
{
  out <- object$spec$data[,c(object$spec$series_id, object$spec$index_id, ".residuals"), with = FALSE]
  return(out)
}

#' Extract Model Sigma
#'
#' @description Extract the standard deviation of the estimated model.
#' @param object an object of class \dQuote{mgam.estimate}.
#' @param ... not currently used.
#' @return A data.table with the index_id, series_id and standard deviation (which
#' can be a vector if the family is \dQuote{gaulss}).
#' @aliases sigma
#' @method sigma mgam.estimate
#' @rdname sigma
#' @export
#'
#'
sigma.mgam.estimate <- function(object, ...)
{
  .residuals <- NULL
  if (object$spec$family$family == "gaulss") {
    out <- object$spec$data[,c(object$spec$series_id, object$spec$index_id, ".sigma"), with = FALSE]
  } else {
    r <- residuals(object)
    s <- r[,list(sd = sd(.residuals)), by = c(object$spec$series_id)]
    out <- data.table(s = object$spec$series, .sigma = as.numeric(s$sd))
    colnames(out)[1] <- object$spec$series_id
  }
  return(out)
}

#' Model Prediction
#'
#' @description Prediction function for class \dQuote{gam.estimate}.
#' @param object an object of class \dQuote{gam.estimate}.
#' @param newdata a data.table of new data to use for the prediction. No missing values allowed.
#' @param nsim The number of simulations to use for generating the simulated predictive distribution.
#' @param copula whether to use a copula to create correlated predictive distribution samples (not available
#' for pooled models). This requires that the newdata are fully aligned in terms of time indices i.e.
#' that the series share the same time indices.
#' @param copula_family the copula family to use.
#' @param shape the degrees of freedom in the case of the Student copula. An NA value means that
#' this will be estimated, but for large dimensional series it may be best to fix this since
#' estimation may be slow.
#' @param distribution a valid distribution from the tsdistribution package which will be used
#' to post-estimate the residuals (marginal) distribution for use in simulating the predictive
#' distribution. In the case of the copula, the quantile function of this distribution will be
#' used to transform the unit hypercube values to innovations.
#' @param tabular whether to whether to return the data in long format as a data.table instead
#' of a list of objects of class tsmodels.predict.
#' @return a list with the following slots:
#' \itemize{
#' \item (tabular = FALSE) prediction_list is a list of \dQuote{tsmodel.predict} objects, one for each series,
#' for which a plot method is available.
#' \item (tabular = TRUE) prediction_table is a long format data.table with columns for by series_id,
#' estimation_date, draw (simulation draw), forecast_date and value (the forecast value).
#' }
#' @note
#' The newdata series_id must match exactly those that were used in the original
#' specification. The function will test for this and optionally reorder the newdata
#' so that the series_id order is the same as in the original dataset. This restriction
#' may be relaxed in a later release so that a pooled model can be used on entirely new
#' series.
#' @aliases predict
#' @method predict mgam.estimate
#' @rdname predict
#' @export
#'
#'
predict.mgam.estimate <- function(object, newdata = NULL, nsim = 5000, copula = TRUE,
                                  copula_family = c("Normal", "Student"), shape = NA,
                                  distribution = NULL, tabular = FALSE,  ...)
{
  # ToDo: return fitted distributional parameters (should also apply to univariate model)
  .std_residuals <- .residuals <- .sigma <- NULL

  # check the newdata id's. Validates that all series_ids are the same as in the model
  # and orders the newdata so that it matches the series_id order of the model
  newdata <- validate_prediction_ids(newdata, object$spec$data, object$spec$series_id)
  copula_family <- match.arg(copula_family[1], c("Normal", "Student"))
  pooled <- object$spec$pooled
  if (copula & pooled) {
    stop("\ncannot use a copula in a pooled model")
  }
  if (is.null(newdata)) {
    stop("\nnewdata cannot be NULL")
  } else {
    if (!is.data.table(newdata)) stop("\nnewdata must be a data.table object")
    if (any(is.na(newdata))) stop("\nNA's found in newdata.")
  }
  check <- validate_prediction_inputs(newdata, object$spec$formula_vars)
  if (!check) stop("\nvariables names of newdata do not match those in the input data")
  check <- validate_prediction_inputs(newdata, object$spec$series_id)
  if (!check) stop(paste0("\nnewdata does not contain the series_id :", object$spec$series_id))
  check <- validate_prediction_inputs(newdata, object$spec$index_id)
  if (!check) stop(paste0("\nnewdata does not contain the index_id :", object$spec$index_id))

  # check that newdata has time and series id
  if (copula) {
    check <- validate_prediction_indices(newdata, object$spec$index_id, object$spec$series_id)
    if (!check$validate) stop('\nnon-matching indices found in newdata between different series. The series should have matching time indices when copula = TRUE')
    h <- as.numeric(check$h)[1]
    .x <- copy(object$spec$data)
    .x[,.std_residuals := .residuals/.sigma]
    residuals_matrix <- dcast(.x, paste0(object$spec$index_id,"~",object$spec$series_id), value.var = ".std_residuals", fun.aggregate = mean)
    residuals_matrix <- residuals_matrix[,-1]
    # scale to z scores
    C <- cor(residuals_matrix, method = "spearman", use = "complete.obs")
    P <- P2p(C)
    if (copula_family == "Normal") {
      cop <- normalCopula(P, dim = ncol(C), dispstr = "un")
    } else if (copula_family == "Student") {
      cop_spec <- tCopula(dim = ncol(C), dispstr = "un", df = NA_real_, df.fixed = FALSE)
      cop <- fitCopula(cop_spec, data = pobs(as.matrix(residuals_matrix)), method = "itau.mpl", start = c(P2p(C), 5),
                       lower = c(rep(-1, ncol(C)), 4.01), upper = c(rep(1, ncol(C), 40)))
    }
    U <- rCopula(nsim * h, cop)
    pset <- split(newdata, by = c(object$spec$series_id))
    gam_models <- tsconvert(object)
    b <- future_lapply(1:length(object$model), function(i){
      innov <- matrix(U[,i], ncol = h, nrow = nsim)
      ndata <- pset[[i]]
      setcolorder(ndata, c(object$spec$index_id, setdiff(names(ndata), object$spec$index_id)))
      ndata <- as.xts(ndata)
      p <- predict(gam_models[[i]], newdata = ndata, nsim = nsim, tabular = FALSE, distribution = distribution, innov = innov, innov_type = "q")
      return(p)
    }, future.packages = c("xts","tsmethods","tsgam","mgcv","data.table"), future.seed = TRUE)
    b <- eval(b)
  } else {
    if (pooled) {
      b <- pooled_prediction(object, newdata, nsim = nsim, tabular = FALSE, distribution = distribution)
    } else {
      b <- nonpooled_prediction(object, newdata, nsim = nsim, tabular = FALSE, distribution = distribution)
    }
  }
  if (tabular) {
    out <- lapply(1:length(b), function(i){
      tabular_predictions(b[[i]]$distribution, object$spec, id = object$spec$series[i])
    })
    out <- list(prediction_table = rbindlist(out))
    xspec <- object$spec
    xspec$data <- NULL
    out$spec <- xspec
    return(out)
  } else {
    xspec <- object$spec
    xspec$data <- NULL
    out <- list(prediction_list = b, spec = xspec)
    return(out)
  }
}

#' Convert an mgam.estimate model to a list of gam.estimate models
#'
#' @description Converts an mgam.estimate model to a list of gam.estimate models
#' @param object an object of class \dQuote{mgam.estimate}.
#' @param ... not currently used.
#' @return A list with each slot of class \dQuote{gam.estimate} holding
#' the estimated model for each series.
#' @aliases tsconvert
#' @method tsconvert mgam.estimate
#' @rdname tsconvert
#' @export
#'
#'
tsconvert.mgam.estimate <- function(object, ...)
{
  if (object$spec$pooled) {
    stop("\nno translation for to gam.estimate for pooled model.")
  }
  eset <- split(object$spec$data, by = c(object$spec$series_id))
  model_list <- lapply(1:length(object$model), function(i) {
    extra <- list()
    extra$time_index <-  eset[[i]][[object$spec$index_id]]
    extra$clean_time_index <-  eset[[i]][[object$spec$index_id]]
    extra$time_class <- object$spec$time_class
    extra$time_zone <- object$spec$time_zone
    out <- list(model = object$model[[i]], spec = extra)
    class(out) <- "gam.estimate"
    return(out)
  })
  return(model_list)
}

pooled_prediction <- function(object, newdata, nsim = 9000, tabular = FALSE, distribution = NULL, ...)
{
  # setup data.table variable names to avoid notes in checking
  prediction <- parameter <- forecast_date <- estimation_date <- draw <- value <- NULL
  `.` <- list
  .fitted <- .sigma <- .residuals <- .mu <- .sigma <- NULL
  newx <- as.data.frame(newdata)
  if (any(is.na(newdata))) stop("\nno NAs allowed in data")
  prediction <- predict(object$model, newdata = newx, type = "response")
  if (is.matrix(prediction)) {
    mean_prediction <- prediction[,1]
    sigma_prediction <- 1/prediction[,2]
    sigma_model <- 1/fitted(object$model ,type = "response")[,2]
  } else {
    mean_prediction <- prediction
    sigma_model <- sd(residuals(object$model, type = "response"))
    sigma_prediction <- rep(sigma_model, length(mean_prediction))
  }
  # create the decomposition
  res <- residuals(object$model, type = "response")
  z <- res/sigma_model
  prediction_table <- copy(newdata)
  prediction_table[,.mu := mean_prediction]
  prediction_table[,.sigma := sigma_prediction]
  prediction_table <- split(prediction_table, by = c(object$spec$series_id))
  actuals_table <- split(object$spec$data, by = c(object$spec$series_id))
  n <- length(object$spec$series)
  if (!is.null(distribution)) {
    if (object$model$family$family != "gaussian" & object$model$family$family != "gaulss") stop("\ndistribution option only available for models estimated using gaussian or gaulss family")
    # estimate distribution on residuals of model
    # predict given mu and distributional parameters
    spec <- distribution_modelspec(z, distribution = distribution)
    distribution_fit <- estimate(spec)
    p_matrix <- distribution_fit$spec$parmatrix
    skew <- p_matrix[parameter == "skew"]$value
    shape <- p_matrix[parameter == "shape"]$value
    lambda <- p_matrix[parameter == "lambda"]$value
    b <- future_lapply(1:n, function(i){
      tmp <- prediction_table[[i]]
      fut_dates <- tmp[[object$spec$index_id]]
      simulated_draws <- do.call(cbind, lapply(1:NROW(tmp), function(j) {
        rdist(distribution = distribution, n = nsim, mu = tmp$.mu[j], sigma = tmp$.sigma[j],
              skew = skew, shape = shape, lambda = lambda)}))
      simulated_draws <- as.matrix(simulated_draws)
      colnames(simulated_draws) <- as.character(fut_dates)
      class(simulated_draws) <- "tsmodel.distribution"
      attr(simulated_draws, "date_class") <- object$spec$time_class
      simulated_draws <- list(original_series = xts(actuals_table[[i]][[object$spec$response_var]], actuals_table[[i]][[object$spec$index_id]]),
                              distribution = simulated_draws,
                              mean = xts(tmp$.mu, fut_dates))
      class(simulated_draws) <- c("gam.predict","tsmodel.predict")
      return(simulated_draws)
    }, future.packages = c("xts","tsmethods","tsgam","mgcv","data.table"), future.seed = TRUE)
    b <- eval(b)
  } else {
    b <- future_lapply(1:n, function(i){
      tmp <- prediction_table[[i]]
      fut_dates <- tmp[[object$spec$index_id]]
      simulated_draws <- as.data.table(predicted_samples(object$model, n = nsim, data = tmp))
      simulated_draws <- dcast(simulated_draws, draw~row, value.var = "response")
      simulated_draws <- simulated_draws[order(draw)]
      simulated_draws <- simulated_draws[,draw := NULL]
      setcolorder(simulated_draws, paste0(1:ncol(simulated_draws)))
      simulated_draws <- as.matrix(simulated_draws)
      colnames(simulated_draws) <- as.character(fut_dates)
      class(simulated_draws) <- "tsmodel.distribution"
      attr(simulated_draws, "date_class") <- object$spec$time_class
      simulated_draws <- list(original_series = xts(actuals_table[[i]][[object$spec$response_var]], actuals_table[[i]][[object$spec$index_id]]),
                              distribution = simulated_draws,
                              mean = xts(tmp$.mu, fut_dates))
      class(simulated_draws) <- c("gam.predict","tsmodel.predict")
      return(simulated_draws)
    }, future.packages = c("xts","tsmethods","tsgam","mgcv","gratia","data.table"), future.seed = TRUE)
    b <- eval(b)
  }
  return(b)
}

nonpooled_prediction <- function(object, newdata, nsim = 9000, tabular = FALSE, distribution = NULL, ...)
{
  # setup data.table variable names to avoid notes in checking
  .mu <- .sigma <- prediction <- parameter <- forecast_date <- estimation_date <- draw <- value <- NULL
  `.` <- list
  if (any(is.na(newdata))) stop("\nno NAs allowed in data")
  prediction <- lapply(1:length(object$model), function(i) {
    newx <- newdata[eval(parse(text = paste0(object$spec$series_id,"==",object$spec$series[i])))]
    out <- predict(object$model[[i]], newdata = newx, type = "response")
    if (is.matrix(out)) {
      mean_prediction <- out[,1]
      sigma_prediction <- 1/out[,2]
      sigma_model <- 1/fitted(object$model[[i]], type = "response")[,2]
    } else {
      mean_prediction <- out
      sigma_model <- sd(residuals(object$model[[i]], type = "response"))
      sigma_prediction <- rep(sigma_model, length(mean_prediction))
    }
    res <- residuals(object$model[[i]], type = "response")
    ret <- list(predicted = data.table(mean = mean_prediction,
                                       sigma = sigma_prediction), actual = data.table(residuals = res, sigma = sigma_model))
    return(ret)
  })
  mean_prediction <- do.call(c, lapply(1:length(prediction), function(i) prediction[[i]]$predicted$mean))
  sigma_prediction <- do.call(c, lapply(1:length(prediction), function(i) prediction[[i]]$predicted$sigma))
  res <- do.call(c, lapply(1:length(prediction), function(i) prediction[[i]]$actual$residuals))
  sigma_model <- do.call(c, lapply(1:length(prediction), function(i) prediction[[i]]$actual$sigma))
  # create the decomposition
  z <- res/sigma_model
  prediction_table <- copy(newdata)
  prediction_table[,.mu := mean_prediction]
  prediction_table[,.sigma := sigma_prediction]
  prediction_table <- split(prediction_table, by = c(object$spec$series_id))
  actuals_table <- split(object$spec$data, by = c(object$spec$series_id))
  n <- length(object$spec$series)
  if (!is.null(distribution)) {
    if (object$model[[1]]$family$family != "gaussian" & object$model[[1]]$family$family != "gaulss") stop("\ndistribution option only available for models estimated using gaussian or gaulss family")
    spec <- distribution_modelspec(z, distribution = distribution)
    distribution_fit <- estimate(spec)
    p_matrix <- distribution_fit$spec$parmatrix
    skew <- p_matrix[parameter == "skew"]$value
    shape <- p_matrix[parameter == "shape"]$value
    lambda <- p_matrix[parameter == "lambda"]$value
    b <- future_lapply(1:n, function(i){
      tmp <- prediction_table[[i]]
      fut_dates <- tmp[[object$spec$index_id]]
      simulated_draws <- do.call(cbind, lapply(1:NROW(tmp), function(j) {
        rdist(distribution = distribution, n = nsim, mu = tmp$.mu[j], sigma = tmp$.sigma[j],
              skew = skew, shape = shape, lambda = lambda)}))
      simulated_draws <- as.matrix(simulated_draws)
      colnames(simulated_draws) <- as.character(fut_dates)
      class(simulated_draws) <- "tsmodel.distribution"
      attr(simulated_draws, "date_class") <- object$spec$time_class
      simulated_draws <- list(original_series = xts(actuals_table[[i]][[object$spec$response_var]], actuals_table[[i]][[object$spec$index_id]]),
                              distribution = simulated_draws,
                              mean = xts(tmp$.mu, fut_dates))
      class(simulated_draws) <- c("gam.predict","tsmodel.predict")
      return(simulated_draws)
    }, future.packages = c("xts","tsmethods","tsgam","mgcv","data.table"), future.seed = TRUE)
    b <- eval(b)
  } else {
    b <- future_lapply(1:n, function(i){
      tmp <- prediction_table[[i]]
      fut_dates <- tmp[[object$spec$index_id]]
      simulated_draws <- as.data.table(predicted_samples(object$model[[i]], n = nsim, data = tmp))
      simulated_draws <- dcast(simulated_draws, draw~row, value.var = "response")
      simulated_draws <- simulated_draws[order(draw)]
      simulated_draws <- simulated_draws[,draw := NULL]
      setcolorder(simulated_draws, paste0(1:ncol(simulated_draws)))
      simulated_draws <- as.matrix(simulated_draws)
      colnames(simulated_draws) <- as.character(fut_dates)
      class(simulated_draws) <- "tsmodel.distribution"
      attr(simulated_draws, "date_class") <- object$spec$time_class
      simulated_draws <- list(original_series = xts(actuals_table[[i]][[object$spec$response_var]], actuals_table[[i]][[object$spec$index_id]]),
                              distribution = simulated_draws,
                              mean = xts(tmp$.mu, fut_dates))
      class(simulated_draws) <- c("gam.predict","tsmodel.predict")
      return(simulated_draws)
    }, future.packages = c("xts","tsmethods","tsgam","mgcv","gratia","data.table"), future.seed = TRUE)
    b <- eval(b)
  }
  return(b)
}

tabular_predictions <- function(x, spec, id)
{
  series_id <- estimation_date <- draw <- NULL
  forecast_date <- value <- NULL
  `.` <- list
  tabular_results <- as.data.table(unclass(x))
  tabular_results[,draw := 1:.N]
  tabular_results <- melt(tabular_results, id.vars = "draw", measure.vars = 1:(ncol(tabular_results) - 1), variable.name = "forecast_date", value.name = "value")
  dclass <- attr(x, "date_class")
  if (dclass == "POSIXct") {
    fun <- fastPOSIXct
  } else {
    fun <- fastDate
  }
  tabular_results[,forecast_date := fun(forecast_date, tz = spec$time_zone)]
  e_date <- max(spec$data[eval(parse(text = paste0(spec$series_id,"==", id)))][[spec$index_id]], na.rm = TRUE)
  tabular_results[,estimation_date := e_date]
  tabular_results[,series_id := id]
  tabular_results <- tabular_results[,list(series_id, estimation_date, draw, forecast_date, value)]
  return(tabular_results)
}


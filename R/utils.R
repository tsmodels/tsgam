validate_prediction_inputs <- function(newdata, variables)
{
  newdata_names <- colnames(newdata)
  check <- all(variables %in% newdata_names)
  return(check)
}

validate_prediction_ids <- function(newdata, olddata, series_id)
{
  newdata_id <- sort(unique(newdata[[series_id]]))
  olddata_id <- sort(unique(olddata[[series_id]]))
  if (!all.equal(newdata_id, olddata_id)) stop("\nnewdata series_id do not match those in the model")
  olddata_id <- unique(olddata[[series_id]])
  xnewdata <- lapply(1:length(olddata_id), function(i){
    return(newdata[eval(parse(text = paste0(series_id,"==",olddata_id[i])))])
  })
  xnewdata <- rbindlist(xnewdata)
  return(xnewdata)
}

validate_prediction_indices <- function(newdata, index_id, series_id)
{
  .indicator <- NULL
  test <- copy(newdata)
  test[,.indicator := 1]
  tab <- dcast(test, paste0(index_id, "~", series_id), value.var = ".indicator", fun.aggregate = length, fill = FALSE, drop = FALSE)
  if (any(is.na(tab))) {
    check <- list(validate = FALSE, h = NULL)
    return(check)
  } else {
    cnames <- colnames(tab)[-1]
    D <- tab[, lapply(.SD, "max"), .SDcols = cnames]
    h <- tab[, lapply(.SD, "sum"), .SDcols = cnames]
    if (any(D > 1)) message("\nvalidation check on newdata indices found duplicate time indices.")
    check <- list(validate = TRUE, h = h)
    return(check)
  }
}

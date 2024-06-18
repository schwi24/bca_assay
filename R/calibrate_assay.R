#' Calibration curves from standard samples
#'
#' This function is used for absolute quantification using a calibration curve 
#' and takes in the standard samples. It extracts the fitted model, the model's 
#' parameters, and Efron's coefficients of determination (pseudo R-squared)
#' including the adjusted version for multiple parameters.
#' @param standards A data frame with the standard samples. It must contain 
#' columns named `A560`, and `Quantity` with numerical values.
#' @param low_cutoff A numerical value or vector with the lower cut-off boundary
#' for the `Quantity` that is used for the calibration. The default is `0`.
#' @param high_cutoff A numerical value or vector with the upper cut-off 
#' boundary for the `Quantity` that is used for the calibration. The default is 
#' `Inf`.
#' @param interval A string with the type of interval to be calculated. Can be 
#' one of "confidence" (default), "prediction" or "none".
#' @param level A numeric value between 0 and 1 that indicates the confidence 
#' level for the calculated intervals (if any). The default is 0.95 for 
#' confidence intervals.
#' 
#' @return A list of tibbles named `fit`, `parameters`, and `r_squared`, 
#' followed by the `model` object.
#' 
#' The tibble in `fit` can be used for plotting the calibration curve with 
#' confidence intervals. It contains the columns `quantity`, `fit`, `upr`, and 
#' `lwr`.
#' The tibble in `parameters` contains the columns `parameter`, `estimate`, 
#' `std_error`, `t_value`, and `pr_t` for each model parameter from the summary
#' of the fitted model.
#' The tibble `r_squared` contains the columns `pseudo_r_squared` and 
#' `adj_r_squared` with Efron's pseudo R-squared and its adjusted version for 
#' multiple fitted parameters.
#' 
#' @importFrom dplyr as_tibble
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom investr predFit
#' @importFrom janitor clean_names
#' @importFrom minpack.lm nlsLM
#' @importFrom tibble add_column
#' @importFrom tibble tibble
#' @export
fit_logistic_calibration_curve <- function(
  standards,
  low_cutoff = 0,
  high_cutoff = Inf,
  interval = "confidence",
  level = 0.95
){
  
  ## Apply cutoffs
  data <- janitor::clean_names(dat = standards) |>
    dplyr::filter(
      .data = _,
      quantity >= low_cutoff,
      quantity <= high_cutoff
    )
    
  
  ## Data range for simulate values
  mn <- min(dplyr::pull(data, quantity))
  mx <- max(dplyr::pull(data, quantity))
  conc <- tibble::tibble(quantity = seq(from = mn, to = mx, length.out = 100))
  
  ## Starting values
  a0 <- dplyr::filter(
    .data = data,
    quantity == min(dplyr::pull(.data = data, quantity))
  ) |>
    dplyr::pull(.data = _, a560)
  d0 <- dplyr::filter(
    .data = data,
    quantity == max(dplyr::pull(.data = data, quantity))
  ) |>
    dplyr::pull(.data = _, a560)
  b0 <- 1
  c0 <- (d0 - a0) / 2
  
  ## Model
  model <- try(
    minpack.lm::nlsLM(
      formula = a560 ~ ((a - d)/(1 + (quantity/c)^b) + d),
      data = data,
      start = list(a = a0, b = b0, c = c0, d = d0),
      lower = c(0, -Inf, 0, 0),
      upper = c(Inf, Inf, Inf, Inf)
    )
  )
  
  ## Interval
  interval <- try(
    investr::predFit(
      object = model,
      newdata = conc,
      interval = "confidence",
      se.fit = FALSE,
      level = level
    )
  )
  
  ## Parameters
  mod_summary <- try(
    summary(model)
  )
  
  ## Pseudo R squared
  pseudo_r2 <- try(
    .pseudo_r_squared(model)
  )
  
  
  ## Results
  mod <- NULL
  if(!inherits(model, "try-error")){
    mod <- model
  }
  fit <- NULL
  if(!inherits(interval, "try-error")){
    fit <- dplyr::as_tibble(interval) |>
      tibble::add_column(.data = _, quantity = conc$quantity, .before = 1)
  }
  params <- NULL
  if(!inherits(mod_summary, "try-error")){
    params <- mod_summary$parameters |>
      janitor::clean_names(dat = _) |>
      dplyr::as_tibble(x = _, rownames = "parameter")
  }
  r_squared <- NULL
  if(!inherits(pseudo_r2, "try-error")){
    r_squared <- pseudo_r2
  }
  
  res <- list(
    fit = fit,
    parameters = params,
    r_squared = r_squared,
    model = mod
  )
  
  return(res)
  
}


#' Quantify Unknowns with Calibration Curve
#'
#' @param data a tibble with columns `quantity` and `a560`
#' @param calibration_curve a calibration curve as created by 
#' `fit_logistic_calibration_curve`
#'
#' @return a tibble with updated values in quantity based on the calibration 
#' curve
#' @export
#' @importFrom janitor clean_names
#' @importFrom tidyr pivot_wider
quantify_unknowns <- function(data, calibration_curve){
  
  params <- calibration_curve$parameters[, c("parameter", "estimate")] |>
    pivot_wider(data = _, names_from = "parameter", values_from = "estimate") |>
    as.list(x = _)
  a <- params$a
  b <- params$b
  c <- params$c
  d <- params$d
  
  data <- janitor::clean_names(dat = data)
  unkowns <- !(data[,"task"] == "Standard" | data[,"task"] == "Empty")
  data[unkowns, "quantity"] <- c * 
    ((a - d)/(data[unkowns, "a560"] - d) - 1)^(1/b)
    
  return(data)

}



#' Efron's R-squared
#'
#' This function calculates Efron's coefficients of determination including an
#' adjusted version for multiple parameters.
#' Based on https://stackoverflow.com/questions/14530770/calculating-r2-for-a-nonlinear-least-squares-fit
#' 
#' @param model A object with a linear or non-linear model
#' @return A tibble with pseudo R-squared value and its adjusted version.
#' @importFrom tibble tibble
.pseudo_r_squared <- function(model){
  
  ## https://stackoverflow.com/questions/14530770/calculating-r2-for-a-nonlinear-least-squares-fit

  pred <- predict(model)  # yhat
  
  ## residual sum of squares (rss)
  n <- length(pred)
  res <- resid(model)  # residuals/errors
  w <- weights(model)
  if (is.null(w)) w <- rep(1, n)
  rss <- sum(w * res ^ 2)
  
  ## total sum of squares (tss)
  resp <- pred + res  # response = yhat + error/residuals
  center <- weighted.mean(resp, w)
  tss <- sum(w * (resp - center)^2)
  
  ## pseudo r2 and adjusted pseudo r2
  r_sq <- 1 - rss/tss
  r_df <- summary(model)$df[2]  # degree of freedoms
  int_df <- 1
  adj_r_sq <- 1 - (1 - r_sq) * (n - int_df) / r_df
  
  ## results
  res <- tibble::tibble(
    pseudo_R_squared = r_sq,
    adj_R_squared = adj_r_sq
  )
  return(res)
  
}


#' Select smoothing parameter for smoothing splines based on leave-one-out CV error
#'
#' @param x predictor variable.
#' @param y response variable.
#' @param cv_lambda vector of candidate lambda values, must be between 0 and 1.
#'
#' @return lamdba value that minimizes leave-one-out CV error.
#' @export
sel_smoothing_para <- function(x, y, cv_lambda) {
  cv_size = length(cv_lambda)
  cv_error <- array(0, dim=cv_size)
  for (i in 1:cv_size) {
    solution <- css_train(x, y, cv_lambda[i])
    cv_error[i] <- cal_loo_cv_error(y, css_predict(x_test = x, basis = solution), solution$S)
  }

  df <- as.data.frame(x = cv_lambda)
  df[,2] <- cv_error
  best_lambda <- df$cv_lambda[which.min(apply(df$V2,MARGIN=1,min))]  # get best lambda

  return(best_lambda)
  #return(df)
}


#' Calculte leave-one-out CV error for a linear smoother
#'
#' @param y response variable values
#' @param f_hat fitted response variable values
#' @param S smoother matrix
#'
#' @return leave-one-out cross-validation error
cal_loo_cv_error <- function(y, f_hat, S) {
  residual <- y - f_hat
  normalizing_weights <- 1 - diag(S)
  cv_error <- mean((residual/normalizing_weights) ^ 2)

  return(cv_error)
}
#' Select smoothing parameter for smoothing splines based on leave-one-out CV error
#'
#' @param x predictor variable.
#' @param y response variable.
#' @param cv_lambda vector of candidate lambda values, must be between 0 and 1.
#'
#' @return lamdba value that minimizes leave-one-out CV error.
#' @export
#' @examples
#' set.seed(1)
#' x_train <- seq(1, 10, 0.1)
#' y_train <- cos(x_train)^3 * 3 - sin(x_train)^2 * 2 + x_train + exp(1)+rnorm(length(x_train),0,1)
#' plot(x_train,y_train)
#' x_test <- seq(1, 10, 0.1)
#' lines(x_test,cos(x_test)^3 * 3 - sin(x_test)^2 * 2 + x_test + exp(1),col="red")
#'
#' cv_lambda <- seq(0,0.001,0.0001)
#' results <- sel_smoothing_para(x_train, y_train, cv_lambda)
#' plot(results$df$cv_lambda, results$df$error, type="o")
#' title('LOO CV Error of Different Lambdas')
#' abline(v = results$best, col='red', lty=2)
#' legends <- c("LOO CV Error", "Best Lambda")
#' legend('topleft', legend=legends, col=c('black', 'red'), lty=c(1,2), cex=0.8)
sel_smoothing_para <- function(x, y, cv_lambda) {
  cv_size = length(cv_lambda)
  cv_error <- array(0, dim=cv_size)
  for (i in 1:cv_size) {
    solution <- css_train(x, y, cv_lambda[i])
    cv_error[i] <- cal_loo_cv_error(y, css_predict(x_test = x, basis = solution), solution$S)
  }

  df <- as.data.frame(x = cv_lambda)
  df['error'] <- cv_error
  best_lambda <- df$cv_lambda[which.min(apply(df$error,MARGIN=1,min))]  # get best lambda

  #return(best_lambda)
  return(list('df' = df, 'best' = best_lambda))
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

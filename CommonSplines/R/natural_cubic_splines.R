#' Regression using natural cubic splines
#'
#' This function provides regression using natural cubic splines with truncated power basis functions.
#' Only univariate input can be used.
#'
#' @param x_train The input vector of training dataset.
#' @param y_train The output vector of training dataset.
#' @param x_test The input values at which evaluations are required.
#' @param df Degrees of freedom. One can supply df rather than knots;
#' ncs() then chooses (df + 1) knots at uniform quantiles of x.
#' The default, df = 4, sets 5 knots with 3 inner knots at uniform quantiles of x.
#' @param knots Breakpoints that define the spline.
#' The default is five knots at uniform quantiles c(0, .25, .5, .75, 1).
#' Typical values are the mean or median for one knot, quantiles for more knots.
#'
#' @return
#' \item{y_pred}{A vector of dimension length(x), the prediction vector evaluated at x_test values.}
#'
#' @export
#'
#' @examples
#' x_train <- seq(1, 10, 0.1)
#' y_train <- cos(x_train)^3 * 3 - sin(x_train)^2 * 2 + x_train + exp(1)+rnorm(length(x_train),0,1)
#' plot(x_train,y_train)
#' title('Comparison of Different Degrees of Freedom')
#' x_test <- seq(1, 10, 0.1)
#' lines(x_test,cos(x_train)^3 * 3 - sin(x_train)^2 * 2 + x_train + exp(1),col="red")
#'
#' df <- 2
#' y_pred <- ncs(x_train, y_train, x_test, df)
#' lines(x_test,y_pred, col='blue')
#' df <- 4
#' y_pred <- ncs(x_train, y_train, x_test, df)
#' lines(x_test,y_pred, col='green')
#' df <- 10
#' y_pred <- ncs(x_train, y_train, x_test, df)
#' lines(x_test,y_pred, col='black')
#' legends <- c("Actual", "Prediction: 2 df", "Prediction: 4 df", "Prediction: 10 df")
#' legend('topleft', legend=legends, col=c('red', 'blue', 'green', 'black'), lty=1, cex=0.8)
ncs <- function(x_train, y_train, x_test, df = NULL, knots = NULL)
{
  train_result <- ncs_train(x_train, y_train,
                                              df = df, knots = knots)
  y_pred <- ncs_predict(x_test, train_result$betas,
                                          train_result$knots, train_result$nknots)
  return(y_pred)
}


#' Generate an evaluated basis matrix for natural cubic splines
#'
#' @param x_train The input vector of training dataset.
#' @param y_train The output vector of training dataset.
#' @param df Degrees of freedom. One can supply df rather than knots;
#' ncs() then chooses (df + 1) knots at uniform quantiles of x.
#' The default, df = 4, sets 5 knots with 3 inner knots at uniform quantiles of x.
#' @param knots Breakpoints that define the spline, in terms of quantiles of x.
#' The default is five knots at uniform quantiles c(0, .25, .5, .75, 1).
#' Typical values are the mean or median for one knot, quantiles for more knots.
#'
#' @return A list of following components:
#' \item{nknots}{Number of knots.}
#' \item{knots}{A vector of knot locations.}
#' \item{N}{Basis matrix evaluated at each x value.}
#' \item{betas}{Least sqaure fit parameters.}
#' @export
#'
#' @examples
#' x_train <- seq(1, 10, 0.1)
#' y_train <- cos(x_train)^3 * 3 - sin(x_train)^2 * 2 + x_train + exp(1)+rnorm(length(x_train),0,1)
#' plot(x_train,y_train)
#' x_test <- seq(1, 10, 0.1)
#' df <- 10
#' train_result <- ncs_train(x_train, y_train, df)
#' print(train_result$betas)
#' print(train_result$N[1:5,1:5])
ncs_train <- function(x_train, y_train, df = NULL, knots = NULL)
{
  ginv <- MASS::ginv
  # get all necessary spline properties
  if (is.null(df) & is.null(knots)) {  # neither is specified
    df <- 4  # Default 4 degrees of freedom
    nknots <- df + 1
    knots <- place_knots(nknots, x_train)
  } else if (is.null(df)) {  # knots is specified
    nknots <- length(knots)
    knots <- quantile(x_train, knots, type=1)
  } else if (is.null(knots)) {
    nknots <- df + 1
    knots <- place_knots(nknots, x_train)
  }

  # evaluate basis functions and obtain basis matrix
  N <- ncs_eval_basis(x_train, knots, nknots)

  # least sqaure fit
  betas <- (ginv(t(N)%*%N))%*%t(N)%*%y_train

  return(list('nknots' = nknots, 'knots' = knots,
              'N' = N, 'betas' = betas))
}


#' Prediction based on trained regression model
#'
#' @param x_test The input values at which evaluations are required.
#' @param betas Least sqaure fit parameters obtained from training.
#' @param knots Knots location in terms of quantiles of x_train, optional, default will be evenly
#' spaced quantiles based on number of knots.
#' @param nknots Number of knots used in training.
#'
#' @return
#' \item{y_pred}{A vector of dimension length(x), the prediction vector evaluated at x_test values.}
#' @export
ncs_predict <- function(x_test, betas, knots, nknots)
{
  N_test = ncs_eval_basis(x_test, knots, nknots)
  y_pred = N_test %*% betas

  return(y_pred)
}


#' Find evenly spaced knots by quantile
#'
#' Knots found include boundary knots at 0th and 100th quantile.
#'
#' @param nknots Number of knots to be located.
#' @param x Data vector on which knots are placed.
#'
#' @return A named vector with knot quantiles and values.
#' @export
place_knots <- function(nknots, x)
{
  if(nknots > 0L) {
    knots <- seq.int(from = 0, to = 1,
                     #length.out = nknots + 2L)[-c(1L, nknots + 2L)]
                     length.out = nknots)
    knots <- quantile(x, knots, type=1)  # type=1 for using inverse of empirical cdf
  }
  return(knots)
}


#' Evaluate basis functions as each x and return the evaluated basis matrix N
#'
#' @param x Predictor variable vector.
#' @param knots Knots location in terms of quantiles of x_train, optional, default will be evenly
#' spaced quantiles based on number of knots.
#' @param nknots Number of knots useded in training.
#'
#' @return Basis matrix evaluated at each x value.
ncs_eval_basis <- function(x, knots, nknots)
{
  N <- matrix(0, nrow=length(x), ncol=nknots)
  for (m in 1:length(x)){
    for (n in 1:nknots){
      # evaluate each basis function
      N[m, n] <- ncs_eval_single_basis(x[m], n, knots, nknots)
    }
  }

  return(N)
}


# Evalute x based on truncated power basis functions for natural cubic splines
#
# @param x A single predictor variable value
# @param i Location index for x vector.
# @param knots Knot location vector.
# @param nknots Number of knots useded in training.
#
# @return Basis function evaluation at x
ncs_eval_single_basis <- function(x, i, knots, nknots)
{
  if (i == 1) {
    return(1)
  } else if (i == 2) {
    return(x)
  } else {
    return(d_k_function(x, k=i-2, knots, nknots)
           - d_k_function(x, k=nknots-1, knots, nknots))
  }
}


# k = 1, ..., K-2, where K = nknots
d_k_function <- function(x, k, knots, nknots)
{
  if (k == nknots) {
    return(0)
  } else {
    d_k = (max(0,x-knots[k])^3 - max(0,x-knots[nknots])^3) / (knots[nknots]-knots[k])
    return(d_k)
  }
}


#' Select smoothing parameter based on leave-one-out CV error
#'
#' @param x predictor variable
#' @param y response variable
#' @param cv_lambda vector of candidate lambda values
#'
#' @return lamdba value that minimizes leave-one-out CV error
#' @export
sel_smoothing_para <- function(x, y, cv_lambda) {
  cv_size = length(cv_lambda)
  cv_error <- array(0, dim=cv_size)
  for (i in 1:cv.size) {
    data <- smoothingSplineTrain(x, y, cv_lambda[i])
    cv_error[i] <- cal_loo_cv_error(y, data$f_hat, data$S)
  }

  df <- as.data.frame(x = cv_lambda)
  df[,2] <- cv_error
  best_lambda <- df$cv_lambda[which.min(apply(df,MARGIN=1,min))]  # get best lambda

  return(best_lambda)
}


#' Calculte leave-one-out CV error
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






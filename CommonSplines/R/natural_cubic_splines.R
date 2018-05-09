#' Regression using natural cubic splines
#'
#' This function provides regressions using natural cubic splines with truncated power basis functions.
#' Only univariate input can be used.
#'
#' @param x_train The input vector of training dataset.
#' @param y_train The output vector of training dataset.
#' @param x_test The input values at which evaluations are required.
#' @param df The degree of freedom specified by user, number of knots will be equal to df.
#' @param knots Knots location in terms of quantiles of x_train, optional, default will be evenly
#' spaced quantiles based on number of knots
#'
#' @return
#' \item{y_pred}{A vector of dimension length(x)The prediction vector evaluated at x_test values}
#'
#' @export
#'
#' @examples
#' x_train <-seq(0, 1, 0.001)
#' y_train <- x^3 * 3 - x^2 * 2 + x + exp(1)+rnorm(length(x),0,0.1)
#' plot(x,y)
#' df <- 10
#' x_test <- seq(0, 1, 0.01)
#' y_pred <- natural_cubic_splines(x, y, x_test, df)
#' plot(x_test,y_pred)
#' lines(x_test,x_test^3 * 3 - x_test^2 * 2 + x_test + exp(1),col="red")
natural_cubic_splines <- function(x_train, y_train, x_test, df = NULL, knots = NULL)
{
  train_result <- natural_cubic_splines.train(x_train, y_train, df = df)
  y_pred <- natural_cubic_splines.predict(x_test, train_result$betas,
                                          train_result$knots, nknots = df)
  return(y_pred)
}


#' Generate an evaluated basis matrix for natural cubic splines
#'
#' @param x_train The input vector of training dataset.
#' @param y_train The output vector of training dataset.
#' @param df The degree of freedom specified by user, number of knots will be equal to df.
#' @param knots Knots location in terms of quantiles of x_train, optional, default will be evenly
#' spaced quantiles based on number of knots
#' @param intercept Default false, do not change.
#'
#' @return A list of following components:
#' \item{knots} A vector of knot locations
#' \item{N} Basis matrix evaluated at each x value
#' \item{betas} Least sqaure fit parameters
#' @export
#'
#' @examples
#' x_train <-seq(0, 1, 0.001)
#' y_train <- x^3 * 3 - x^2 * 2 + x + exp(1)+rnorm(length(x),0,0.1)
#' plot(x,y)
#' df <- 10
#' x_test <- seq(0, 1, 0.01)
#' train_result <- natural_cubic_splines.train(x, y, df)
#' print(train_result$betas)
#' print(train_result$N[1:5,1:5])
natural_cubic_splines.train <- function(x_train, y_train, df = NULL, knots = NULL,
                                        intercept = FALSE)
{
  # get all necessary spline properties
  nknots <- df
  knots <- place_knots(nknots, x_train)

  # evaluate basis functions and obtain basis matrix
  N <- eval_basis_functions(x_train, knots, nknots)

  # least sqaure fit
  betas <- (solve(t(N)%*%N))%*%t(N)%*%y_train

  return(list('knots' = knots, 'N' = N, 'betas' = betas))
}


#' Prediction based on trained regression model
#'
#' @param x_test The input values at which evaluations are required.
#' @param betas Least sqaure fit parameters obtained from training.
#' @param knots Knots location in terms of quantiles of x_train, optional, default will be evenly
#' spaced quantiles based on number of knots
#' @param nknots Number of knots used in training.
#'
#' @return
#' \item{y_pred}{A vector of dimension length(x)The prediction vector evaluated at x_test values}
#' @export
natural_cubic_splines.predict <- function(x_test, betas, knots, nknots){
  N_test = eval_basis_functions(x_test, knots, nknots)
  y_pred = N_test %*% betas

  return(y_pred)
}


#' Find evenly spaced out knots by quantile
#'
#' @param nknots Number of knots to be located.
#' @param x Data vector on which knots are placed.
#'
#' @return A named vector with knot quantiles and values
#' @export
place_knots <- function(nknots, x) {
  if(nknots > 0L) {
    knots <- seq.int(from = 0, to = 1,
                     length.out = nknots + 2L)[-c(1L, nknots + 2L)]
    quantile(x, knots, type=1)  # type=1 for using inverse of empirical cdf
  }
  return(knots)
}


#' Evaluate basis functions as each x and return the evaluated basis matrix N
#'
#' @param x Predictor variable vector
#' @param knots Knots location in terms of quantiles of x_train, optional, default will be evenly
#' spaced quantiles based on number of knots
#' @param nknots Number of knots useded in training.
#'
#' @return Basis matrix evaluated at each x value
eval_basis_functions <- function(x, knots, nknots) {
  N <- matrix(0, nrow=length(x), ncol=nknots)
  for (m in 1:length(x)){
    for (n in 1:nknots){
      N[m, n] <- basis_function(x[m], n, knots, nknots)  # evaluate each basis function
    }
  }

  return(N)
}


#' Evalute x based on truncated power basis functions for natural cubic splines
#'
#' @param x A single predictor variable value
#' @param i Location index for x vector.
#' @param knots Knot location vector.
#' @param nknots Number of knots useded in training.
#'
#' @return Basis function evaluation at x
basis_function <- function(x, i, knots, nknots)
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
d_k_function <- function(x, k, knots, nknots){
  if (k == nknots) {
    return(0)
  } else {
    d_k = (max(0,x-knots[k])^3 - max(0,x-knots[nknots])^3) / (knots[nknots]-knots[k])
    return(d_k)
  }
}

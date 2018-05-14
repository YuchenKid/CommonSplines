#' Train regression coefficients for natural cubic splines.
#'
#' In the least square fitting of nonparametric regression coefficients,
#' Moore-Penrose generalized inverse (ginv\{MASS\}) is used to aviod computational problems.
#'
#' @param x_train The input vector of training dataset.
#' @param y_train The output vector of training dataset.
#' @param df Degrees of freedom. One can supply df rather than knots;
#' \code{ncs_train} then chooses (df + 1) knots at uniform quantiles of x.
#' The default, df = 4, sets 5 knots with 3 inner knots at uniform quantiles of x.
#' @param knots Breakpoints that define the spline, in terms of quantiles of x or real values of x.
#' The default is five knots at uniform quantiles c(0, .25, .5, .75, 1).
#' Typical values are the mean or median for one knot, quantiles for more knots.
#' @param q A boolean variable indicating whether \code{knots} provided are quantiles or real values. When \code{q}=TRUE, \code{knots}
#' provided are quantiles of x. When \code{q}=FALSE, \code{knots} provided are real values of x. Default is FALSE.
#' @return A list of following components:
#' \item{nknots}{Number of knots.}
#' \item{knots}{A vector of knot locations.}
#' \item{N}{Basis matrix evaluated at each x value.}
#' \item{betas}{Least sqaure fit parameters.}
#' @export
#' @references Friedman, J., Hastie, T., & Tibshirani, R. (2001). The Elements of Statistical Learning (Vol. 1, pp. 337-387). New York: Springer series in statistics.
#' Chapter 5.2.1.
#'
#' Venables, W. N. and Ripley, B. D. (1999) Modern Applied Statistics with S-PLUS. Third Edition. Springer. p.100.
#'
#' @seealso \code{ncs_predict}, \code{np_reg}
#' @examples
#' x_train <- seq(1, 10, 0.1)
#' y_train <- cos(x_train)^3 * 3 - sin(x_train)^2 * 2 + x_train + exp(1)+rnorm(length(x_train),0,1)
#' plot(x_train,y_train)
#' x_test <- seq(1, 10, 0.1)
#' df <- 10
#' train_result <- ncs_train(x_train, y_train, df)
#' print(train_result$beta)
#' print(train_result$N[1:5,1:5])
ncs_train <- function(x_train, y_train, df = NULL, knots = NULL, q = FALSE)
{
  ginv <- MASS::ginv

  knots <- generate_knots(x_train = x_train, df = df, knots = knots, q = q)

  # evaluate basis functions and obtain basis matrix
  N <- ncs_basis(x_train, knots)

  # least sqaure fit
  betas <- (ginv(t(N)%*%N))%*%t(N)%*%y_train

  return(list('knots' = knots,
              'N' = N, 'beta' = betas,'nknots'=length(knots)))
}


#' Prediction using regression by natural cubic splines.
#'
#' @param x_test The input values at which evaluations are required.
#' @param basis The return value of function \code{ncs_train}.
#' Instead of specify \code{knots} and \code{beta},One can supply \code{basis} directly.
#' @param knots Breakpoints that define the spline. \code{knots} should be in terms of real-values of x
#'  It can be the return value of \code{generate_knots}.
#' @param beta The coefficients of nonparametric regression.
#'
#' @return
#' \item{y_pred}{A vector of dimension length(x), the prediction vector evaluated at x_test values.}
#' @seealso \code{ncs_train}, \code{np_reg}
#' @export
ncs_predict <- function(x_test,knots=NULL,beta=NULL,basis=NULL)
{
  if (!is.null(basis)) #basis is not null
  {
    knots <- basis$knots
    beta <- basis$beta
  }
  if (is.null(basis) & (is.null(knots) | is.null(beta))) #basis is null and some parameters are null
  {
    print("Some necessary parameters (knots/beta) are missing!")
  }

  N_test <- ncs_basis(x_test, knots)
  y_pred <- N_test %*% beta

  return(y_pred)
}

#' Generate an evaluated basis matrix for natural cubic splines
#'
#' @param x Predictor variable vector.
#' @param knots Knots location in terms of real values of x.
#'
#' @return Basis matrix evaluated at each x value.
#' @examples
#' x<-seq(0, 1, 0.001)
#' knots <- seq(0, 1, 0.1)
#'
#' basis<-ncs_basis(x,knots)
#' plot(x,rep(0,length(x)),type="l",ylim=c(0,1))
#' for (i in 1: (length(knots))){
#'   lines(x,basis[,i],col=i)
#' }
#' @export
ncs_basis <- function(x, knots)
{
  N <- matrix(0, nrow=length(x), ncol=length(knots))
  for (m in 1:length(x)){
    for (n in 1:length(knots)){
      # evaluate each basis function
      N[m, n] <- ncs_eval_single_basis(x[m], n, knots, length(knots))
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


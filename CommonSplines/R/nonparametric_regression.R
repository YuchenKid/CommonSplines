#' Nonparametric Regression using spline based methods
#'
#' This function provides regression using natural cubic splines with truncated power basis functions.
#' Only univariate input can be used.
#'
#' @param x_train The input vector of training dataset.
#' @param y_train The output vector of training dataset.
#' @param x_test The input values at which evaluations are required.
#' @param df Degrees of freedom. One can supply df rather than knots;
#' (df + 1) knots are chosen at uniform quantiles of x.
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
#' y_pred <- np_reg(x_train, y_train, x_test,func="ncs", df=df)
#' lines(x_test,y_pred, col='blue')
#' df <- 4
#' y_pred <- np_reg(x_train, y_train, x_test,func="ncs", df=df)
#' lines(x_test,y_pred, col='green')
#' df <- 10
#' y_pred <- np_reg(x_train, y_train, x_test,func="ncs", df=df)
#' lines(x_test,y_pred, col='black')
#' legends <- c("Actual", "Prediction: 2 df", "Prediction: 4 df", "Prediction: 10 df")
#' legend('topleft', legend=legends, col=c('red', 'blue', 'green', 'black'), lty=1, cex=0.8)
np_reg<- function(x_train, y_train, x_test, func = 'bs',order=3,df = NULL, knots = NULL,lambda=0.001,q=FALSE)
{
  if(func=="ncs"){
    fit <- ncs_train(x_train, y_train, df = df, knots = knots,q=q)
    y_pred <- ncs_predict(x_test, fit$betas,fit$knots)
  }
  else if(func=="bs"){
    fit <- bs_train(x_train, y_train, real_knots = knots,q=q)
    y_pred <- ncs_predict(x_test,order=fit$order, fit$knots, fit$betas)
  } else if(func=="css"){
    fit <- css_train(x_train, y_train,lambda)
    y_pred <- css_predict(x_test, fit$betas,fit$knots)
  }else if(func=="pbs"){
    fit <- pbs_train(x_train, y_train,order=order, df = df, knots = knots,q=q)
    y_pred <- ncs_predict(x_test, order=fit$order,fit$betas,fit$knots)
  }else{
    print("A correct regression spline/smoothing spline function needs to be specified.")
  }
  return(y_pred)
}
#' Generate knots when real value is not specified.
#'
#' @param x_train The input vector of training dataset.
#' @param df Degrees of freedom. One can supply df rather than knots;
#' \code{generate_knots} then chooses (df + 1) knots at uniform quantiles of x.
#' The default, df = 4, sets 5 knots with 3 inner knots at uniform quantiles of x.
#' @param knots Breakpoints that define the spline, in terms of quantiles or real valus of x.
#' The default is five knots at uniform quantiles c(0, .25, .5, .75, 1).
#' Typical values are the mean or median for one knot, quantiles for more knots.
#' @param q A boolean variable define whether \code{knots} provided are quantiles or real values. When \code{q}=TRUE, \code{knots}
#' provided are quantiles of x. When \code{q}=FALSE, \code{knots} provided are real values of x.

#' @return A vector of knots in terms of real values of x.
#' @export
generate_knots<-function(x_train,df,knots,q)
{
  # get all necessary spline properties
  if (is.null(df) & is.null(knots)) {  # neither is specified
    df <- 4  # Default 4 degrees of freedom
    nknots <- df + 1
    knots <- place_knots(nknots, x_train)
  } else if (is.null(df)&q==TRUE) {  # knots is specified as quantiles of x
    nknots <- length(knots)
    knots <- quantile(x_train, knots, type=1)
  } else if (is.null(df)&q==FALSE) {  # knots is specified as real values of x
    knots<-knots
  } else if (is.null(knots)) {
    nknots <- df + 1
    knots <- place_knots(nknots, x_train)
  }else{
    print("degree of freedom/knots is required!")
  }
  return(knots)
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

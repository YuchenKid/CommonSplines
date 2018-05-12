#' Train a smoothing spline with squared 2nd derivative penalty using natural cubic spline
#'
#' This function trains a smoothing spline with squared 2nd derivative penalty. It has an explicit,finite-dimensional,
#' unique minimizer which is a natural cubic spline. This function can be used for small or moderate number of knots. When the number of data N<=50, all knots are included.
#' When N>50, 50 knots are uniformly chosen from the training dataset.
#'
#' @param x The input vector of training dataset.
#' @param y The output vector of training dataset.
#' @param lambda A fixed smoothing parameter.
#' @return A list with the following components:
#' \item{beta}{ The coefficients of natural splines.}
#' \item{S}{The smoother matrix.}
#' \item{knots}{The knots used to construct the B-splines, including innerknots, boundary knots and phantom knots}
#' @examples
#' x<-seq(0, 1, 0.001)
#' y <- x^3 * 3 - x^2 * 2 + x + exp(1)+rnorm(length(x),0,0.1)
#' plot(x,y)
#' lambda<-0.001
#'
#' basis<-css_train(x,y,lambda)
#' cat("the knots chosen are: ",basis$knots)
#' @export

css_train<-function (x,y,lambda)
{
  basis<- ncs_train(x, y, knots=x/max(x))
  knots<-unname(basis$knots)
  h<-numeric(length(knots)-1)
  for (i in 1:length(h)){
    h[i]<-knots[i+1]-knots[i]
  }
  delta<-matrix(0,nrow=length(knots)-2,ncol=length(knots))
  w<-matrix(0,nrow=length(knots)-2,ncol=length(knots)-2)
  for(i in 1:(length(knots)-2)){
    delta[i,i]<-1/h[i]
    delta[i,i+1]<--1/h[i]-1/h[i+1]
    delta[i,i+2]<-1/h[i+1]
    w[i-1,i]<-h[i]/6
    w[i,i-1]<-h[i]/6
    w[i,i]<-(h[i]+h[i+1])/3
  }
  K<-t(delta)%*%MASS::ginv(w)%*%delta
  I<-diag(length(knots))
  S<-solve(I+lambda*K)%*%y
  beta<-MASS::ginv(basis$N)%*%S
  solution<-list("beta"=beta,"S" = S, "knots"=x)
  return(solution)
}
#' Prediction using smoothing spline with squared 2nd derivative penalty
#'
#' This function takes the coefficients trained by \code{CubicSmoothingSpline.Train} and evaluate the output at x_test
#'
#' @param x_test The input values at which evaluations are required.
#' @param basis The return value of function \code{css_train}.
#' Instead of specify \code{knots} and \code{beta},One can supply \code{basis} directly.
#' @param knots Breakpoints that define the spline. \code{knots} should be in terms of real-values of x
#'  It can be the return value of \code{generate_knots}.
#' @param beta The coefficients of nonparametric regression.
#'
#' @return The evaluated output at x_test.
#' @examples
#' x<-seq(0, 1, 0.0015)
#' y <- x^3 * 3 - x^2 * 2 + x + exp(1)+rnorm(length(x),0,0.1)
#' plot(x,y)
#' lambda<-0.001
#' basis<-css_train(x,y,lambda)
#'
#' x_test<-seq(0, 1, 0.1)
#' fit<-css_predict(x_test=x_test,basis=basis)
#'
#' plot(x_test,fit)
#' lines(x_test,x_test^3 * 3 - x_test^2 * 2 + x_test + exp(1),col="red")

#' @export
css_predict<-function (x_test,knots=NULL,beta=NULL,basis=NULL)
{
  if(is.null(basis)==0) #basis is not null
  {
    knots<-basis$knots
    beta<-basis$beta
  }
  if(is.null(basis)&(is.null(knots)|is.null(beta))) #basis is null and some parameters are null
  {
    print("some necessary parameters (knots/beta) are missing!")
  }
  N_test<-ncs_basis(x_test, knots)
  f<-N_test%*%beta
  return (f)
}

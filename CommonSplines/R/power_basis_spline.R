#' Evaluate basis functions as each x and return the evaluated basis matrix N
#'
#' @param x Predictor variable vector.
#' @param knots Knots location in terms of quantiles of x_train, optional, default will be evenly
#' spaced quantiles based on number of knots.
#' @param nknots Number of knots useded in training.
#'
#' @return Basis matrix evaluated at each x value.
pbs_basis <- function(x,order,knots)
{
  phi <- matrix (0, nrow=length(x), ncol=(length(knots)+order))
  for(m in 1:length(x)){
    for(n in 1:(length(knots)+order)){
      if(n <= order){
        phi[m,n] <- x[m]^(n-1)
      }
      else{
        temp <- (x[m]-knots[n-order])^(order-1)
        phi[m,n] <- ifelse(temp>0,temp,0)
      }
    }
  }
  return(phi)
}

#' Regression using Power Basis spline
#'
#' This function provides regressions using Power Basis splines. The basis are defined as
#' 1,x,x^2,...,x^m,(x-k1)^(m-1)+,(x-k2)^(m-1)+,...,(x-kn)^(m-1)+
#' where m is the order, k1, k2 and kn are n knots, '+' denotes the positive part.
#'
#' Only univariate input can be used.
#'
#' @param x The input vector of training dataset.
#' @param y The output vector of training dataset.
#' @param x_test The input values at which evaluations are required.
#' @param order The order that defines the spline.
#' @param knots The internal knots that define the spline.
#' @return A list with the following components:
#' \item{beta}{ The coefficients of nonparametric regression.}
#' \item{basis}{The spline basis matrix of dimension c(length(x), length(knots)+order)}
#' \item{f}{The evaluated output at x_test.}
#' @examples
#' n <- 100
#' t <- seq(0,2*pi,length.out = 100)
#' a <- 3
#' b <- 2
#' c.unif <- runif(n)
#' amp <- 2
#' set.seed(1)
#' y1 <- a*sin(b*t)+c.unif*amp # uniform error
#' knots <- c(min(t),2*pi*c(1/4,2/4,3/4),max(t))
#' order <- 4
#' basis <- pbs_train(t,y1,order,knots)
#' fit<-pbs_predict(t,basis=basis)
#' y.hat <- fit
#' plot(t, y1, t="l")
#' lines(t, y.hat, col=2)
#' @export
pbs_train <- function(x,y,order,knots)
{
  ginv <- MASS::ginv
  knots <- unique(sort(knots))
  G<-pbs_basis(x,order,knots)
  y <- matrix(y)
  Beta <- (ginv(t(G)%*%G))%*%t(G)%*%y



  solution<-list("beta"=Beta,"basismatrix"=G,"knots"=knots, "order"=order)
  return (solution)
}
pbs_predict<-function (x_test,order=NULL,knots=NULL,beta=NULL,basis=NULL) #knots should contain phantom knots or generated from bs_knots
{
  if(is.null(basis)==0) #basis is not null
  {
    knots<-basis$knots
    order<-basis$order
    beta<-basis$beta
  }
  if(is.null(basis)&(is.null(order)|is.null(knots)|is.null(beta))) #basis is null and some parameters are null
  {
    print("some necessary parameters (order/knots/beta) are missing!")
  }
  phi<-pbs_basis(x_test,order,knots)
  f<-phi%*%beta
  return (f)
}

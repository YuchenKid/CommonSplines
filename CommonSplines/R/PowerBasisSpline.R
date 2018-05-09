#' Regression using Power Basis spline
#'
#' This function is a generalization of CubicPowerBasisSpline with arbitrary order
#' 
#' Only univariate input can be used.
#' 
#' @param x The input vector of training dataset.
#' @param y The output vector of training dataset.
#' @param x_test The input values at which evaluations are required.
#' @param order The order that defines the spline.
#' @param innerknots The internal knots that define the spline.
#' @return A list with the following components:
#' \item{beta}{ The coefficients of nonparametric regression.}
#' \item{basis}{The spline basis matrix of dimension c(length(x), NumKnots+order)}
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
#' innerknots <- 2*pi*c(1/4,2/4,3/4)
#' order <- 4
#' solution <- PowerBasisSpline(t,y1,t,order,innerknots)
#' y.hat <- solution$f
#' plot(t, y1, t="l")
#' lines(t, y.hat, col=2)
#' @export
PowerBasisSpline <- function(x,y,x_test,order,innerknots)
{
  innerknots <- unique(sort(innerknots))
  NumX <- length(x)
  NumKnots <- length(innerknots)
  G <- matrix(0, nrow=NumX, ncol=(NumKnots+order)) # matrix G is NumX*(NumKnots+order) dimensional, used for coefficient estimation.
  y <- matrix(y)
  
  for(m in 1:NumX){
    for(n in 1:(NumKnots+order)){
      if(n <= order){
        G[m,n] <- x[m]^(n-1)
      }
      else{
        temp <- (x[m]-innerknots[n-order])^(order-1)
        G[m,n] <- ifelse(temp>0,temp,0)
      }
    }
  }
  Beta <- (ginv(t(G)%*%G))%*%t(G)%*%y
  
  NumTest <- length(x_test)
  phi <- matrix (0, nrow=NumTest, ncol=(NumKnots+order))
  for(m in 1:NumTest){
    for(n in 1:(NumKnots+order)){
      if(n <= order){
        phi[m,n] <- x_test[m]^(n-1)
      }
      else{
        temp <- (x_test[m]-innerknots[n-order])^(order-1)
        phi[m,n] <- ifelse(temp>0,temp,0)
      }
    }
  }
  f <- phi%*%Beta
  
  solution<-list("beta"=Beta,"basis"=G,"f"=f)
  return (solution)
}

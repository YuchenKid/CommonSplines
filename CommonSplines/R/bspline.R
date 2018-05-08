#' Regression using B-spline basis
#'
#' This function provides nonparametric regressions using B-splines. The B-splines are defined following
#' the recursive formulas due to de Boor. Only univariate input can be used.
#'
#' @param x The input vector of training dataset.
#' @param y The output vector of training dataset.
#' @param x_test The input values at which evaluations are required.
#' @param order The order of B-spline functions. The default is order=4 for cubic B-splines.
#' @param innerknots The internal knots that define the spline.
#' @return A list with the following components:
#' \item{beta}{ The coefficients of nonparametric regression.}
#' \item{basis}{The B-spline basis matrix of dimension c(length(x), df). df = length(innerknots) + order.}
#' \item{f}{The evaluated output at x_test.}
#' @examples
#' x<-seq(0, 1, 0.001)
#' y <- x^3 * 3 - x^2 * 2 + x + exp(1)+rnorm(length(x),0,0.1)
#' plot(x,y)
#'
#' innerknots <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
#' order<-4
#' x_test<-seq(0, 1, 0.01)
#'
#' b_fit<-bspline(x,y,x_test,order,innerknots)
#'

#' plot(x_test,b_fit$f)
#' lines(x_test,x_test^3 * 3 - x_test^2 * 2 + x_test + exp(1),col="red")
#'

#' plot(x,rep(0,length(x)),type="l",ylim=c(0,1))
#' for (i in 1: (j+order)){
#' lines(x,b_fit$basis[,i])
#' }
#' @export
bsplineBasis <-function (x,y,x_test,order=4,innerknots) #degree<-order-1
{
  highknot = max(x,innerknots)
  lowknot = min(x,innerknots)
  innerknots <- unique (sort (innerknots))
  #knots <-c(rep(lowknot, order), innerknots, rep(highknot, order))
  knots <-c(rep(lowknot-0.0001, order-1),lowknot, innerknots,highknot, rep(highknot+0.0001, order-1))
  n <- length (x)
  j <- length (innerknots) # j+order basis functions
  G <- matrix (0,  nrow=n,ncol=(j+2*order-1)) # matrix G is n*(J+q) dimensional, used for coefficient estimation.
  for (i in 1: (j+2*order-1)){
    G[,i]<-ifelse((x>=knots[i]&x<knots[i+1]),1,0)
  }
  for (k in 2:order){
    N<-G
    for (i in 1: (j+2*order-k)){
      if((knots[i+(k-1)]-knots[i])==0&&(knots[i+k]-knots[i+1])==0){
        G[,i]<-rep(0,n)
      }else if ((knots[i+(k-1)]-knots[i])==0){
        G[,i]<-(knots[i+k]-x)/(knots[i+k]-knots[i+1])*N[,(i+1)]
      }else if ((knots[i+k]-knots[i+1])==0){
        G[,i]<-(x-knots[i])/(knots[i+(k-1)]-knots[i])*N[,i]
      }else{
        G[,i]<-(x-knots[i])/(knots[i+(k-1)]-knots[i])*N[,i]+(knots[i+k]-x)/(knots[i+k]-knots[i+1])*N[,(i+1)]
      }
    }
  }
  Gb<-G[,1:(j+order)]
  beta<-(solve(t(Gb)%*%Gb))%*%t(Gb)%*%y

  m<-length(x_test)
  phi <- matrix (0,  nrow=m,ncol=(j+2*order-1))
  for (i in 1: (j+2*order-1)){
    phi[,i]<-ifelse((x_test>=knots[i]&x_test<knots[i+1]),1,0)
  }
  for (k in 2:order){
    N<-phi
    for (i in 1: (j+2*order-k)){
      if((knots[i+(k-1)]-knots[i])==0&&(knots[i+k]-knots[i+1])==0){
        phi[,i]<-rep(0,m)
      }else if ((knots[i+(k-1)]-knots[i])==0){
        phi[,i]<-(knots[i+k]-x_test)/(knots[i+k]-knots[i+1])*N[,(i+1)]
      }else if ((knots[i+k]-knots[i+1])==0){
        phi[,i]<-(x_test-knots[i])/(knots[i+(k-1)]-knots[i])*N[,i]
      }else{
        phi[,i]<-(x_test-knots[i])/(knots[i+(k-1)]-knots[i])*N[,i]+(knots[i+k]-x_test)/(knots[i+k]-knots[i+1])*N[,(i+1)]
      }
    }
  }
  phi_b<-phi[,1:(j+order)]
  f<-phi_b%*%beta
  solution<-list("beta"=beta ,"basis" = Gb,"f"=f)
  return (solution)
}


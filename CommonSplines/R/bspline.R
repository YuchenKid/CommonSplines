#' Add phantom knots for B-splines
#'
#' @param x Predictor variable vector.
#' @param real_knots The innerknots and boundary knots that define the spline. The knots can all be innerknots.
#' The knots provided can be quantiles of x or real values.
#' More explanation of \code{knots}, \code{df}, \code{q} can be seen in \code{generate_knots}.
#' @param df Degrees of freedom. One can supply df rather than knots.
#' @param q A boolean variable define whether knots provided are quantiles or real values. When \code{q}=TRUE, \code{real_knots}
#' are quantiles of x. When \code{q}=FALSE, \code{real_knots} are real values of x. Default is FALSE.
#' @param order The order of basis functions. order=degree+1
#' @return The knots used to construct the B-splines, including innerknots, boundary knots and phantom knots.
#' @export
bs_knots<-function(x,df=NULL,real_knots=NULL,q=FALSE,order)
{
  knots<-generate_knots(x,df,real_knots,q)
  highknot = max(x)
  lowknot = min(x)
  knots <- unique (sort (knots))
  knots <-c(rep(lowknot-0.0001, order), knots, rep(highknot+0.0001, order))
  return(knots)
}
#' Generate an evaluated basis matrix for B-splines
#'
#' This function generates B-spline basis. The B-splines are defined following
#' the recursive formulas due to de Boor. Only univariate input can be used.
#'
#' @param x Predictor variable vector.
#' @param knots The knots used to construct the B-splines, including innerknots, boundary knots and phantom knots.
#' It can be generated by \code{bs_knots}.
#' @param order The order of basis functions. order=degree+1
#' @return Basis matrix evaluated at each x value.
#' @references De Boor, C., De Boor, C., Mathématicien, E. U., De Boor, C., & De Boor, C. (1978). A Practical Guide to Splines (Vol. 27, p. 325). New York: Springer-Verlag.
#' @examples
#' x<-seq(0, 1, 0.001)
#' knots <- seq(0, 1, 0.1)
#' order<-4
#' knots<-bs_knots(x,real_knots=knots,order=order)
#'
#' basis<-bs_basis(x,order,knots)
#' plot(x,rep(0,length(x)),type="l",ylim=c(0,1))
#' for (i in 1: (length(knots)-order)){
#'   lines(x,basis[,i],col=i)
#' }
#' @export
bs_basis <-function (x,order,knots)
{
  j <- length (knots)-2*order # j+order basis functions
  n <- length (x)
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
  return (G[,1:(j+order)])
}


#' Train regression coefficients for B-splines.
#'
#' @param x The input vector of training dataset.
#' @param y The output vector of training dataset.
#' @param order The order of B-spline functions. The default is order=4 for cubic B-splines.
#' @param real_knots The innerknots and boundary knots that define the spline.
#' Phantom knots should not be included. Phantom knots will be generated by \code{bs_knots}
#' The knots provided can be quantiles of x or real values.
#' More explanation of \code{knots}, \code{df}, \code{q} can be seen in function \code{generate_knots}.
#' @param df Degrees of freedom. One can supply df rather than knots.
#' @param q A boolean variable define whether knots provided are quantiles or real values. When \code{q}=TRUE, \code{real_knots}
#' are quantiles of x. When \code{q}=FALSE, \code{real_knots} are real values of x. Default is FALSE.
#' @return A list with the following components:
#' \item{beta}{The coefficients of nonparametric regression.}
#' \item{basis}{The B-spline basis matrix of dimension c(length(x), df). df = length(innerknots) + order.}
#' \item{knots}{The knots used to construct the B-splines, including innerknots, boundary knots and phantom knots}
#' \item{order}{The order of basis functions. order=degree+1}
#' @seealso \code{bs_knots}, \code{bs_basis}, \code{bs_predict},\code{np_reg}
#' @references Friedman, J., Hastie, T., & Tibshirani, R. (2001). The Elements of Statistical Learning (Vol. 1, pp. 337-387). New York: Springer series in statistics.
#' Chapter 5, Appendix.
#' @examples
#' x<-seq(0, 1, 0.001)
#' y <- x^3 * 3 - x^2 * 2 + x + exp(1)+rnorm(length(x),0,0.1)
#' plot(x,y)
#' knots <- seq(0, 1, 0.1)
#' order<-4
#'
#' basis<-bs_train(x,y,order,knots)
#' plot(x,rep(0,length(x)),type="l",ylim=c(0,1))
#' for (i in 1: (length(knots)+order)){
#'   lines(x,basis$basismatrix[,i],col=i)
#' }
#' @export
bs_train <-function (x,y,order,real_knots=NULL,df = NULL,q=FALSE) #degree<-order-1
{
  knots<-bs_knots(x=x,df=df,real_knots=real_knots,q=q,order=order)
  G<-bs_basis(x,order,knots)
  beta<-(MASS::ginv(t(G)%*%G))%*%t(G)%*%y
  solution<-list("beta"=beta ,"basismatrix" = G,"knots"=knots, "order"=order)
  return (solution)
}
#' Prediction using regression spline with B-spline basis
#'
#' This function provides prediction at value of interest using regression spline with B-spline basis.
#' The B-splines are generated by \code{bs_basis} and trained by the \code{bs_train}.
#' The return value of \code{bs_train} can be used as an argument of \code{bs_predict}
#'
#' @param x_test The input values at which evaluations are required.
#' @param basis The return value of function \code{bs_train}.
#' Instead of specify \code{knots}, \code{order} and \code{beta},One can supply \code{basis} directly.
#' @param knots Breakpoints that define the spline. \code{knots} should be in terms of real-values of x
#'  and contain innner, boundary and phantom knots. It can be the return value of \code{bs_knots}.
#' @param order The order of basis functions. order=degree+1
#' @param beta The coefficients of nonparametric regression.
#' @return The evaluated output at x_test.
#' @seealso \code{bs_basis},\code{bs_knots}, \code{bs_train}, \code{np_reg}.
#' @examples
#' x<-seq(0, 1, 0.001)
#' y <- x^3 * 3 - x^2 * 2 + x + exp(1)+rnorm(length(x),0,0.1)
#' plot(x,y)
#' knots <- seq(0.1, 0.9, 0.01)
#' order<-4

#' basis<-bs_train(x,y,order,knots)
#'
#' x_test<-seq(0, 1, 0.01)
#' fit<-bs_predict(x_test,basis=basis)

#' plot(x_test,fit)
#' lines(x_test,x_test^3 * 3 - x_test^2 * 2 + x_test + exp(1),col="red")
#' @export
bs_predict<-function (x_test,order=NULL,knots=NULL,beta=NULL,basis=NULL) #knots should contain phantom knots or generated from bs_knots
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
  phi<-bs_basis(x_test,order,knots)
  f<-phi%*%beta
  return (f)
}



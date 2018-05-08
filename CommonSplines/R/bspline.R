#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
bspline <-function (x,y,x_test,order,innerknots) #degree<-order-1
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


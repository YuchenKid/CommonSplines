% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power_basis_spline.R
\name{pbs_basis}
\alias{pbs_basis}
\title{Evaluate basis functions at each x and return the evaluated basis matrix N}
\usage{
pbs_basis(x, order, knots)
}
\arguments{
\item{x}{The input vector of training dataset.}

\item{order}{The order that defines the truncated power basis spline.}

\item{knots}{The innerknots and boundary knots that define the spline. The knots should be real values.
The knots can be generated by \code{generate_knots}.}
}
\value{
Basis matrix evaluated at each x value.
}
\description{
Evaluate basis functions at each x and return the evaluated basis matrix N
}
\examples{
x<-seq(0, 1, 0.001)
knots <- seq(0, 1, 0.1)
order<-4

basis<-pbs_basis(x,order,knots)
plot(x,rep(0,length(x)),type="l",ylim=c(0,1))
for (i in 1: (length(knots)+order)){
  lines(x,basis[,i],col=i)
}
}
\seealso{
\code{generate_knots}.
}

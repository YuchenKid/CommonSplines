#' Generate an evaluated basis matrix for natural cubic splines
#'
#' @param x_train
#' @param y_train
#' @param df
#' @param knots
#' @param intercept
#' @param Boundary.knots
#'
#' @return an evaluated basis matrix of size (length(x) * df)
#' @export
#'
#' @examples
natural_cubic_splines.train <- function(x_train, y_train, df = NULL, knots = NULL,
                                        intercept = FALSE, Boundary.knots = range(x_train))
{
  # get all necessary spline properties
  nknots <- df
  knots <-
    if(nknots > 0L) {
      knots <- seq.int(from = 0, to = 1,
                       length.out = nknots + 2L)[-c(1L, nknots + 2L)]
      quantile(x_train, knots, type=1)  # type=1 for using inverse of empirical cdf
    }

  N <- matrix(0, nrow=length(x_train), ncol=nknots) # basis matrix

  # evaluate basis functions as each x
  for (m in 1:length(x_train)){
    for (n in 1:nknots){
      # evaluate each basis_function(x, i)
      N[m, n] <- basis_function(x_train[m], n)
    }
  }

  # least sqaure fit
  betas <- (solve(t(N)%*%N))%*%t(N)%*%y_train

  return(list('N' = N, 'betas' = betas))
}

#' Evalute x based on truncated power basis functions for natural cubic splines
#'
#' @param x
#' @param i,i = 1, ..., length(x)
#'
#' @return
basis_function <- function(x, i)
{
  if (i == 1){
   return(1)
  } else if (i == 2){
   return(x)
  } else {
    return(d_k_function(x, k=i-2) - d_k_function(x, k=nknots-1))
  }
}

# k = 1, ..., K-2, where K = nknots
d_k_function <- function(x, k){
  if (k == nknots){
    return(0)
  } else {
    d_k = (max(0,x-knots[k])^3 - max(0,x-knots[nknots])^3)/(knots[nknots]-knots[k])
    return(d_k)
  }
}


#' Prediction based on trained regression model
#'
#' @param x_test
#'
#' @return
#' @export
#'
#' @examples
natural_cubic_splines.predict <- function(x_test, betas){
  # Evaluate splines at all x_test
  # Return y_predict values




}




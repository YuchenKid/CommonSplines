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
natural_cubic_splines.train <- function(x_train, y_train, df = NULL, knots = NULL, intercept = FALSE,
               Boundary.knots = range(x))
{
  # get all necessary spline properties
  nknots <- df
  knots <-
    if(nknots > 0L) {
      knots <- seq.int(from = 0, to = 1,
                       length.out = nknots + 2L)[-c(1L, nknots + 2L)]
      quantile(x, knots, type=1)  # type=1 for computation by inverse of empirical cdf
    }

  # evaluate basis functions as each x
  for (x_i in x){
    for (i in 1:nknots){
      # evaluate each basis_function(x, i)
      basis_function(x_i, i)
    }
  }

  # construct predictor basis matrix
  # basis_m = Matrix(length(x), nknots)

  # least sqaure fit
  # fit = lm(y_train ~ (basis_m))

}

#' Evalute x based on truncated power basis functions for natural cubic splines
#'
#' @param x
#' @param i
#'
#' @return
basis_function <- function(x, i)
{
  if (i == 1){
   1
  } else if (i == 2){
   x
  } else {
    d_k_function(k=i-2, x)
  }

}

d_k_function <- function(k, x){
  if (k == nknots){
    0
  } else {
    # d_k = ((x-knot(k))+**3 - (x-knot(K)+**3))/(knot(K)-knot(k))
    print('Not yet lah!')
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
natural_cubic_splines.predict <- function(x_test){
  # Evaluate splines at all x_test
  # Return y_predict values
}

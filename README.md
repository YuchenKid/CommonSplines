# CommonSplines
This is a R package that covers commonly seen nonparametric regression using spline-based methods. For regression spline, commonly seen basis functions are provided such as truncated power basis, natural cubic spline basis, and B-spline basis. For regularization, penalties on squared second-order derivative are provided, i.e., cubic smoothing spline. This package mainly refers to:
> Friedman, J., Hastie, T., & Tibshirani, R. (2001). The elements of statistical learning (Vol. 1, pp. 337-387). New York: Springer series in statistics. Chapter 5.

Below are the steps to install it and try it out. 

1. Install the package from RStudio console. You would need to install “devtools” package if you do not have it already installed.
```
devtools::install_github("YuchenKid/CommonSplines/CommonSplines”)
```
2. Load the library.
```
library(CommonSplines)
```
3. Now you can see our package in the “Packages” tab of RStudio. Inside some of our functions' documentation, you can find a few examples to try out:
```
np_reg()  # Main function for regression, user can choose among several spline methods.
pbs_train()  # Train regression coefficients for truncated power basis splines.
ncs_train()  # Train regression coefficients for natural cubic splines.
bs_train()  # Train regression coefficients for B-splines.
css_train()  # Train regression coefficients for cubic smoothing splines.
sel_smoothing_para()  # Select smoothing parameter for smoothing splines based on leave-one-out CV error.
```
There are more functions for you to find out and play with!

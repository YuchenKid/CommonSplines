# CommonSplines
This is a R package that covers commonly seen nonparametric regression using spline-based methods. For regression spline, commonly seen basis functions are provided such as truncated power basis, natural cubic spline basis, and B-spline basis. For regularization, penalties on squared second-order derivative are provided, i.e., cubic smoothing spline. This package mainly refer to “Friedman, J., Hastie, T., & Tibshirani, R. (2001). The elements of statistical learning (Vol. 1, pp. 337-387). New York: Springer series in statistics,” Chapter 5.

Below are the steps to install it and try it out. 

(1) Install the package from RStudio console. You would need to install “devtools” package if you do not have it already installed.
      devtools::install_github("YuchenKid/CommonSplines/CommonSplines”)
(2) Load the library.
      library(CommonSplines)
Now you can see our package in the “Packages” tab of RStudio. Inside some of our functions' documentation, you can find a few examples to try out:
      np_reg()
      pbs_train()
      ncs_train()
      bs_train()
      css_train()
      sel_smoothing_para()
There are more functions for you to find out and play with!


<!-- README.md is generated from README.Rmd. Please edit that file -->

# TRSmoothLogistic

We provide tensor regression and tensor regression with smooth or lasso
regularization.

## Installation

You must install the following packages first before install
`TRSmoothLogistic`:

``` r
install.packages("ADMM")
# install.packages("devtools")
devtools::install_github("gigi881118/ADMM_adjustment")
```

And you can install the `TRSmoothLogistic` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gigi881118/Tensor-smooth-logistic-regression")
```

## Available Functions

We provide following classes of functions. For more details, please see
help pages of each function using `help()` function in your R session.

|    Function    | Description                                                 |
|:--------------:|:------------------------------------------------------------|
| `interpolate`  | Resize array                                                |
|  `KTtotensor`  | Turn the tensor in decomposition form to the regular tensor |
|  `TensorReg`   | Tensor regression                                           |
| `TensorRegNet` | Tensor regression with Smooth or Lasso Regularization       |

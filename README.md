
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `CSMSN`

<!-- badges: start -->

[![R-CMD-check](https://github.com/SergeNIANGORAN/CSMSN/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SergeNIANGORAN/CSMSN/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of CSMSN is to provide tools for the implementation of the
differentes methods of centralized statistical monitoring for a
continuous variable. This consists of comparisons of the variable
distribution between one center to the distribution of the others
centers in order to detect center with atypical pattern.

## Overview

`CSMSN` is an `R` package that provides differents methods for the
implementation of Centralized Statistical Monitoring (CSM) for
continuous variables data in a multicenter clinical trial. This package
provides one function for generating data following the gaussian
distribution in a multicenter clinical trial with a shift in the mean of
one center, and four functions which performs the CSM method of
**Desmet**, **Hatayama**, **Student** and **Distance**.

The main function of this package which compiles and performs the
simulations of all the four methods together is
`MASTER_CSM_MOY_GLOBAL_SIMS()`.

The methods implemented in this package are detailed in the following
article:

<!-- ## Installation -->
<!-- You can install the development version of CSMSN from [GitHub](https://github.com/) with: -->
<!-- ``` r -->
<!-- # install.packages("devtools") -->
<!-- devtools::install_github("SergeNIANGORAN/CSMSN") -->
<!-- ``` -->
<!-- ## Example -->
<!-- This is a basic example which shows you how to solve a common problem: -->
<!-- ```{r example} -->
<!-- library(CSMSN) -->
<!-- ## basic example code -->
<!-- ``` -->

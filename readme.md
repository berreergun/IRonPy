
<!-- README.md is generated from README.Rmd. Please edit that file -->

IRon
====

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

A Python Package for solving **I**mbalanced **R**egressi**on** tasks

Introduction
------------

Imbalanced domain learning has almost exclusively focused on solving classification tasks, where the objective is to predict cases labelled with a rare class accurately. Such a well-defined approach for regression tasks lacked due to two main factors. First, standard regression tasks assume that each value is equally important to the user. Second, standard evaluation metrics focus on assessing the performance of the model on the most common cases. This package contains methods to tackle imbalanced domain learning problems in regression tasks, where the objective is to predict extreme (rare) values.

The methods contained in this package are: - an automatic and non-parametric method to obtain such relevance functions, building on the concept of relevance proposed by Torgo and Ribeiro (2007); - visualisation tools; - suite of evaluation measures for optimisation/validation processes; - the squared-error relevance area measure, an evaluation metric tailored for imbalanced regression tasks.

Requirements
------------
For windows python32 is needed.

Installation
------------

<!-- You can install the last released version from [CRAN](https://CRAN.R-project.org) with: -->
<!--``` {r, eval=FALSE} -->
<!-- install.packages("IRon") -->
<!-- ``` -->
<!-- You can also install the latest development version from GitHub: -->
<!-- ```{r, eval=FALSE} -->
<!-- install.packages("devtools") -->
<!-- devtools::install_github("nunompmoniz/IRon") -->
<!-- ``` -->
To install from github use the following command lines:

``` r
pip install IRon  # You need to install this package!

```



References
----------

-   Rita Ribeiro and Nuno Moniz (2020). "Imbalanced Regression and Extreme Value Prediction". Machine Learning. 109, pp 1803–1835, Springer. https://doi.org/10.1007/s10994-020-05900-9
-   Luís Torgo and Rita Ribeiro (2007). "Utility-based regression". Proceedings of ECML/PKDD 2007, J. N. Kok, J. Koronacki, R. Lopez deMantaras, S. Matwin, D. Mladenic, and A. Skowron, Eds., pp. 597–604

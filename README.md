# BivarIntCensored

The analysis of the joint cumulative distribution function (CDF) with bivariate event time data is a challenging problem. This package provides a user-friendly toolkit that estimates the joint CDF for bivariate interval-censored data (also current status data) through spline-based seive nonparametric maximum likelihood estimators (NPMLE). Besides, a nonparametric correlation test is provided to check whether two outcomes are correlated or not.

> The details of the method can be found in papers [*A spline-based nonparametric analysis for interval-censored bivariate survival data*](http://www3.stat.sinica.edu.tw/ss_newpaper/SS-2019-0296_na.pdf) and [*Partially monotone tensor spline estimation of the joint distribution function with bivariate current status data*](https://projecteuclid.org/journals/annals-of-statistics/volume-40/issue-3/Partially-monotone-tensor-spline-estimation-of-the-joint-distribution-function/10.1214/12-AOS1016.full).

## Installation
To install the package: 
```r
devtools::install_github("junyzhou10/BivarIntCensored")
```

## Usage
Use main function `BiIntCensd(dat, ...)` to conduct analysis. `dat` is a list containing two items with interval censored data. For details please check the helping documents.

There are two sample simulated data in the package, namely `SampleDat_case1` and `SampleDat_case2`, where 'case1' has current status data for both outcomes and 'case2' has interval-censored data for both outcomes. Please feel free to play with the data (e.g., see examples in the helping document). Notably, this package supports the case that one outcome is current status data and the other is interval-censored data as well. 

## Outputs
The outputs include **nonparametric test statistics** for the correlation between two outcomes as well as the **estimated NPMLE** of the I-spline coefficients. The **surface/curves** of joint and marginal CDF will be plotted automatically. If user provides the time points `[T1, T2]` for the prediction of probability of the events, the estimated marginal and joint CDF will be returned, indicating the chance of observing event marginally or jointly at the given time points.

#' Draw from a bimodal normal distribution
#'
#' Simulates random values from a two-component normal (Gaussian) mixture
#' distribution. With probability `prop`, values are drawn from
#' N(mean1, sd1^2); otherwise from N(mean2, sd2^2).
#'
#' This function is adapted from `rbinorm()` in the FamilyRank package.
#'
#' @param n Integer. Number of observations to simulate.
#' @param mean1 Numeric. Mean of component 1.
#' @param mean2 Numeric. Mean of component 2.
#' @param sd1 Non-negative numeric. Standard deviation of component 1.
#' @param sd2 Non-negative numeric. Standard deviation of component 2.
#' @param prop Numeric in [0, 1]. Mixing probability for component 1.
#'
#' @return A numeric vector of length `n`.
#' @importFrom stats rnorm

rbinorm <-
  function(n, mean1, mean2, sd1, sd2, prop){
    if(prop>1 || prop< 0){stop("proportion should be between 0 and 1")}
    if(n<1){stop("n must be greater than or equal to 1")}
    if(sd1 < 0 || sd2 < 0){stop("standard deviations must be non-negative")}
    z <- rbinom(n,size=1,prob=prop)
    z * rnorm(n, mean=mean1, sd=sd1) + (1-z)*rnorm(n, mean=mean2, sd=sd2)
  }

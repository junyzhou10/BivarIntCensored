#' @title Simulate Bivariate Interval-censored Data (case 2)
#' @description This function generates the bivariate event time data from Clayton copula model, where the marginal distributions of each outcome follows an exponential distribution. The observations windows are generated from uniform distribution.
#' @param size Sample size, i.e., the number of subjects
#' @param lambda1 Hazard (rate parameter) of the exponential distribution for outcome 1
#' @param lambda2 Similar to lambda1, but just for outcome 2
#' @param K.tau Kendall's \eqn{\tau} for Clayton copula. It determines the correlation between two outcomes in the joint distribution
#' @param T1.start The observation time (left boundary) can be no earlier than this value for outcome 1
#' @param T1.end The observation time (right boundary) can be no late than this value for outcome 1
#' @param T2.start Similar to T1.start, for outcome 2
#' @param T2.end Similar to T1.end, for outcome 2
#' @details Changing lambda1 and lambda2 can lead to different proportions of left-, interval,
#' and right-censoring. Another way is to adjust T.start/T.end.
#' Increasing T.start can increase the proportion of left-censoring (which will lead to less interval censoring);
#' while decreasing T.end can increase the proportion of right-censoring (and lower that of interval censoring)
#' @return A list of two elements: T1 and T2. Both including tu_, tv_, and Ind_.
#' [tu, tv] is the observation interval with tu can be as small as 0 (most likely indicating a left-censoring)
#' and tv can be Inf indicating a right-censoring. Ind is the indicator of censoring type.
#' @examples
#' SampleDat_case2 = BiCase2_SimuDat(size = 1000, lambda1 = 0.45, lambda2 = 0.55, K.tau = 0.25,
#' T1.start= 0.1, T1.end = 4, T2.start = 0.1, T2.end = 4)
#' @export
#'

BiCase2_SimuDat <- function(size, lambda1, lambda2, K.tau, T1.start, T1.end, T2.start, T2.end){
  #window.generate.method=1: simple way. if else, then use a complex one
  theta = 2*K.tau/(1-K.tau)
  u = runif(size)
  z = runif(size)
  if ( theta != 0) {
    v = ((z^(-theta/(1+theta)) - 1)*u^(-theta) + 1)^(-1/theta)
    T1 = -log(1 - u)/lambda1
    T2 = -log(1 - v)/lambda2
  } else {
    T1 = -log(1 - u)/lambda1
    T2 = -log(1 - z)/lambda2
  }

  ###############   Generate observing time U & V for two events   ##################
  ## First generate observation windows
  pos1    = seq(ceiling(T1.start), round(T1.end)-1, 1)
  pos2    = seq(ceiling(T2.start), round(T2.end)-1, 1)
  delta1  = seq(0.1, 0.3, length.out =  length(pos1))
  delta2  = seq(0.1, 0.3, length.out =  length(pos2))
  mis.pr1 = seq(0.05, 0.2, length.out =  length(pos1))
  mis.pr2 = seq(0.05, 0.2, length.out =  length(pos2))
  # Here for T1
  window1 = rnorm(size*length(pos1), mean = pos1, sd = delta1)
  window1 = pmax(T1.start, window1)
  # we need to make sure each row is monotone increasing sequence
  window1 = matrix(window1, ncol = length(pos1), byrow = T)
  window1 = as.vector(apply(window1, 1, sort))
  # Make sure the last window is within [T1.start, T1.end]
  window1[window1 > T1.end] = NA
  orig.w1 = matrix(window1, ncol = length(pos1), byrow = T)
  mis.po1 = rbinom(size*length(pos1), size = 1, prob = mis.pr1)
  window1[mis.po1==1] = NA
  mismat1 = matrix(window1, ncol = length(pos1), byrow = T)
  # Now for T2
  window2 = rnorm(size*length(pos2), mean = pos2, sd = delta2)
  window2 = pmax(T2.start, window2)
  # we need to make sure each row is monotone increasing sequence
  window2 = matrix(window2, ncol = length(pos2), byrow = T)
  window2 = as.vector(apply(window2, 1, sort))
  # Make sure the last window is within [T2.start, T2.end]
  window2[window2 > T2.end] = NA
  orig.w2 = matrix(window2, ncol = length(pos2), byrow = T)
  mis.po2 = rbinom(size*length(pos2), size = 1, prob = mis.pr2)
  window2[mis.po2==1] = NA
  mismat2 = matrix(window2, ncol = length(pos2), byrow = T)
  # if a row is all na, then directly replace it by orginal data
  mismat1[rowSums(is.na(mismat1)) == length(pos1),] = orig.w1[rowSums(is.na(mismat1)) == length(pos1),]
  mismat2[rowSums(is.na(mismat2)) == length(pos2),] = orig.w2[rowSums(is.na(mismat2)) == length(pos2),]

  ## Put Two dependent time series T1/T2 into the Observation window1/2 to generate tu_1/2, tv_1/2
  # First T1
  temp1 = mismat1 - T1
  temp1[temp1 < 0] = Inf
  temp1[is.na(temp1)] = Inf
  tv_1 = apply(temp1, 1, min) + T1

  temp1 = T1 - mismat1
  temp1[temp1 < 0] = Inf
  temp1[is.na(temp1)] = Inf
  tu_1 = T1 - apply(temp1, 1, min)
  tu_1[!is.finite(tu_1)] = 0
  Ind_1 = ifelse(tu_1 == 0, "Left", ifelse(is.finite(tv_1), "Interval", "Right"))
  UV_1 = data.frame(tu_1, tv_1, Ind_1)

  # Then T2
  temp1 = mismat2 - T2
  temp1[temp1 < 0] = Inf
  temp1[is.na(temp1)] = Inf
  tv_2 = apply(temp1, 1, min) + T2

  temp1 = T2 - mismat2
  temp1[temp1 < 0] = Inf
  temp1[is.na(temp1)] = Inf
  tu_2 = T2 - apply(temp1, 1, min)
  tu_2[!is.finite(tu_2)] = 0
  Ind_2 = ifelse(tu_2 == 0, "Left", ifelse(is.finite(tv_2), "Interval", "Right"))
  UV_2 = data.frame(tu_2, tv_2, Ind_2)

  #id.boot = matrix(sample(size, size*n.boot, replace = T), nrow = n.boot, byrow = T)

  return(list("T1"=UV_1, "T2"= UV_2))
}

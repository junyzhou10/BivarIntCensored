#' @title OptLogLik (not for call)
#' @description Optimize the loglikelihood of interval censored data
#' @param T1 Observed intervals for outcome 1
#' @param T2 Observed intervals for outcome 2
#' @param l Order of I-spline (degree + 1)
#' @param knot Knot positions for the I-spines
#' @param p_n Number of parameters corresponding to the I-spline basis (for outcome 1)
#' @param q_n Number of parameters corresponding to the I-spline basis (for outcome 2)
#' @return A list of outcomes including the
#' \item{rho}{estimated \eqn{\rho} indicating the correlations between two outcomes}
#' \item{seive.M_ij}{Seive MLE for Mu}
#' \item{seive.w_i}{Seive MLE for w}
#' \item{seive.p_j}{Seive MLE for p}
#'

OptLogLik <- function(T1, T2, l, knot, p_n, q_n) {
  # Generate I-spline data
  Is.u1 = suppressWarnings(iSpline(T1$tu_1, knots=knot$knot1, Boundary.knots=knot$boundary1, degree=l-2, intercept=T))
  Is.v1 = suppressWarnings(iSpline(T1$tv_1, knots=knot$knot1, Boundary.knots=knot$boundary1, degree=l-2, intercept=T))
  Is.u2 = suppressWarnings(iSpline(T2$tu_2, knots=knot$knot2, Boundary.knots=knot$boundary2, degree=l-2, intercept=T))
  Is.v2 = suppressWarnings(iSpline(T2$tv_2, knots=knot$knot2, Boundary.knots=knot$boundary2, degree=l-2, intercept=T))

  LogLikBivar <- function(Mu, w, p) {
    # Part 1: both left censoring
    I11 = T1$Ind_1=="Left" & T2$Ind_2=="Left"
    if (sum(I11)==0) {
      M1  = 0
    } else {
      U1   = Is.v1[I11,,drop = FALSE] # notice the reason here is V not U as in the paper because the data structure we used to save interval data
      U2   = Is.v2[I11,,drop = FALSE]
      M1   = sum(log(diag(U1%*%Mu%*%t(U2))))
      # M1   = sum(log(sapply(seq(nrow(U1)), function(x) sum(U1[x,] %*% t(U2[x,]) * Mu)))) # cannot be recognized due to sapply & crossprod
    }

    # Part 2: left & interval
    I12 = T1$Ind_1=="Left" & T2$Ind_2=="Interval"
    if (sum(I12)==0) {
      M2  = 0
    } else {
      U1  = Is.v1[I12,,drop = FALSE]
      U2  = Is.u2[I12,,drop = FALSE]
      V2  = Is.v2[I12,,drop = FALSE]
      M2  = sum(log(diag(U1%*%Mu%*%t(V2-U2))))
    }

    # Part 3: left & right
    I13 = T1$Ind_1=="Left" & T2$Ind_2=="Right"
    if (sum(I13)==0) {
      M3  = 0
    } else {
      U1  = Is.v1[I13,,drop = FALSE]
      V2  = Is.u2[I13,,drop = FALSE]
      M3  = sum(log(diag(U1%*%Mu%*%t(matrix(1, nrow(V2), ncol(V2)) - V2)) + U1%*%w))
    }

    # Part 4: Interval & Left
    I21 = T1$Ind_1=="Interval" & T2$Ind_2=="Left"
    if (sum(I21)==0) {
      M4  = 0
    } else {
      U1  = Is.u1[I21,,drop = FALSE]
      V1  = Is.v1[I21,,drop = FALSE]
      U2  = Is.v2[I21,,drop = FALSE]
      M4  = sum(log(diag((V1-U1)%*%Mu%*%t(U2))))
    }

    # Part 5: Interval & Interval
    I22 = T1$Ind_1=="Interval" & T2$Ind_2=="Interval"
    if (sum(I22)==0) {
      M5  = 0
    } else {
      U1  = Is.u1[I22,,drop = FALSE]
      V1  = Is.v1[I22,,drop = FALSE]
      U2  = Is.u2[I22,,drop = FALSE]
      V2  = Is.v2[I22,,drop = FALSE]
      M5  = sum(log(diag( (V1-U1)%*%Mu%*%t(V2-U2) )))
    }

    # Part 6: Interval & Right
    I23 = T1$Ind_1=="Interval" & T2$Ind_2=="Right"
    if (sum(I23)==0) {
      M6  = 0
    } else {
      U1  = Is.u1[I23,,drop = FALSE]
      V1  = Is.v1[I23,,drop = FALSE]
      V2  = Is.u2[I23,,drop = FALSE]
      M6  = sum(log(diag((V1-U1)%*%Mu%*%t(matrix(1, nrow(V2), ncol(V2)) - V2)) + (V1-U1)%*%w))
    }

    # Part 7: Right & Left
    I31 = T1$Ind_1=="Right" & T2$Ind_2=="Left"
    if (sum(I31)==0) {
      M7  = 0
    } else {
      V1  = Is.u1[I31,,drop = FALSE]
      U2  = Is.v2[I31,,drop = FALSE]
      M7  = sum(log(diag(U2%*%t(Mu)%*%t(matrix(1,nrow(V1),ncol(V1)) - V1)) + U2%*%p))
    }

    # Part 8: Right & Interval
    I32 = T1$Ind_1=="Right" & T2$Ind_2=="Interval"
    if (sum(I32)==0) {
      M8  = 0
    } else {
      V1  = Is.u1[I32,,drop = FALSE]
      U2  = Is.u2[I32,,drop = FALSE]
      V2  = Is.v2[I32,,drop = FALSE]
      M8  = sum(log(diag((V2-U2)%*%t(Mu)%*%t(matrix(1, nrow(V1), ncol(V1)) - V1)) + (V2-U2)%*%p))
    }

    # Part 9: Right & Right
    I33 = T1$Ind_1=="Right" & T2$Ind_2=="Right"
    if (sum(I33)==0) {
      M9  = 0
    } else {
      V1  = Is.u1[I33,,drop = FALSE]
      V2  = Is.u2[I33,,drop = FALSE]
      M9  = sum(log(1 + diag(V1%*%Mu%*%t(V2)) - V2%*%(t(Mu)%*%matrix(1,p_n-1,1) + p) - V1%*%(Mu%*%matrix(1,q_n-1,1) + w) ) )
    }

    return(M1+M2+M3+M4+M5+M6+M7+M8+M9)
  }

  Mu = CVXR::Variable(p_n-1, q_n-1)
  w  = CVXR::Variable(p_n-1)
  p  = CVXR::Variable(q_n-1)
  # define log-likelihood function that need to be optimized
  objective <- Minimize(-LogLikBivar(Mu, w, p))
  Opt.Problem <- Problem(objective, list(Mu >= 0, w>=0, p>=0, sum(Mu)+sum(w)+sum(p) - 1 <= 0 ))
  CVXR_result <- solve(Opt.Problem)
  seive.M_ij = CVXR_result$getValue(Mu); seive.M_ij[which(seive.M_ij<0,arr.ind = T)] = 0
  seive.w_i  = CVXR_result$getValue(w); seive.w_i[seive.w_i<0] = 0
  seive.p_j  = CVXR_result$getValue(p); seive.p_j[seive.p_j<0] = 0

  # Integral interval (to calculate rho)
  alpha  = 0.01
  upper1 = quantile(c(T1$tu_1[T1$tu_1 != 0], T1$tv_1[is.finite(T1$tv_1)]), 1-alpha)
  lower1 = quantile(c(T1$tu_1[T1$tu_1 != 0], T1$tv_1[is.finite(T1$tv_1)]), alpha)
  upper2 = quantile(c(T2$tu_2[T2$tu_2 != 0], T2$tv_2[is.finite(T2$tv_2)]), 1-alpha)
  lower2 = quantile(c(T2$tu_2[T2$tu_2 != 0], T2$tv_2[is.finite(T2$tv_2)]), alpha)
  # Integration & rho.hat:
  int.BS1 = ibs(upper1, knots=knot$knot1, Boundary.knots=knot$boundary1, degree=l-1, intercept=T) -
    ibs(lower1, knots=knot$knot1, Boundary.knots=knot$boundary1, degree=l-1, intercept=T)
  int.BS2 = ibs(upper2, knots=knot$knot2, Boundary.knots=knot$boundary2, degree=l-1, intercept=T) -
    ibs(lower2, knots=knot$knot2, Boundary.knots=knot$boundary2, degree=l-1, intercept=T)

  int.IS1 = rev(cumsum(rev(int.BS1[2:p_n])))
  int.IS2 = rev(cumsum(rev(int.BS2[2:q_n])))
  F12     = t(int.IS1)%*%seive.M_ij%*%t(t(int.IS2))
  F1      = sum((rowSums(seive.M_ij)+seive.w_i)*int.IS1)
  F2      = sum((colSums(seive.M_ij)+seive.p_j)*int.IS2)
  rho.hat = F12 - F1*F2
  return(list(rho.hat=rho.hat, seive.M_ij=seive.M_ij, seive.w_i=seive.w_i, seive.p_j=seive.p_j))
}

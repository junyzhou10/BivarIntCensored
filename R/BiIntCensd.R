#' @title Bivariate Interval-Censored Data Analysis
#' @description This package provide a tookit which estimates the joint distribution function (as well as marginal distributions of each outcome) for interval-censored data using a spline-based sieve nonparametric MLE
#' @param dat Observed interval-censored data which can be either case 1 (current status data) or case 2 (interval censored data) for either outcome. See details.
#' @param l Order of I-spline (degree + 1). The default is \code{l=4}
#' @param pred.times Input times using format [T1, T2] to predict the probability of observing events marginally and jointly, i.e., [F1, F2, F12]. Multiple inputs of times are supported by providing a n x 2 matrix of times
#' @param int.k_1 User specified number of knots for outcome 1. If not specified, \eqn{N^{1/3}} knots will be chosen at quantiles
#' @param int.k_2 Similar to int.k_1, just for outcome 2
#' @param t1_l (outcome 1) lower bound for nonparametric association test. Default is `NULL` which will automatically adopt \eqn{\tau_{1,l}}
#' @param t1_h (outcome 1) upper bound for nonparametric association test. Default is `NULL` which will automatically adopt \eqn{\tau_{1,h}}
#' @param t2_l (outcome 2) lower bound for nonparametric association test. Default is `NULL` which will automatically adopt \eqn{\tau_{2,l}}
#' @param t2_h (outcome 2) upper bound for nonparametric association test. Default is `NULL` which will automatically adopt \eqn{\tau_{2,h}}
#' @param Corr.Test Whether the test of correlation between two outcomes should be conducted. Default is FALSE. If TRUE, bootstrap will be adopted to calculate \eqn{SE(\rho)}.
#' @param nBootstrp Number of bootstrap samples
#' @param seed Random seed for bootstrap
#' @param parallel Parallel computing, default is TRUE. It is designed for speed up the correlation test. Does not have any effect when \code{Corr.Test == FALSE}
#' @details \code{dat} should meet the following format: \cr
#'     A list of two matrix with name 'T1' and 'T2'. Each one is a n x 3 matrix with first column for left boundaries of the observation interval, second column for right boundaries,
#'     and third column for indicators of censoring type, in c("Left", "Right", "Interval"). If current status data, there is no "Interval".\cr\cr
#'     Notice: If there are only two columns provided, the package assumes that they are for interval times [T1, T2] and the indicators will be generated automatically following the rule: \cr \cr
#'     If T1==0 -> "Left" censored \cr
#'     If is.infinite(T2) -> "Right" censored \cr
#'     Else -> "Interval" censored \cr
#' @return A list object including:
#' \item{rho.hat}{The estimated \eqn{\rho}, implies the correlation between two outcomes. See reference papers for more details}
#' \item{MC.rho}{If \code{Corr.Test = TRUE}, it yields a vector of all estimated \eqn{\rho} from bootstraps}
#' \item{seive.M_ij, seive.w_i, seive.p_j}{Estimated I-spline coefficients. See reference papers for more details}
#' \item{F1.hat, F2.hat, F12.hat}{Estimated marginal and joint CDFs. They mainly used for making the plot of the CDFs}
#' \item{Pred.Probs}{If \code{pred.times} is given, it contains the estimated [F1(T1), F2(T2), F12(T1,T2)] at given time points}
#' @author Junyi Zhou, Yuan Wu, and Ying Zhang.
#' @references Wu, Y., & Zhang, Y. (2012). Partially monotone tensor spline estimation of the joint distribution function with bivariate current status data. \emph{The Annals of Statistics, 40(3)}, 1609-1636.\cr \cr
#'     Wu Y., Zhang, Y., & Zhou, J. \emph{Statistica Sinica Preprint No: SS-2019-0296}.
#' @examples
#' res = BiIntCensd(SampleDat_case1,
#'                  Corr.Test = FALSE,
#'                  pred.times = rbind(c(1,2), c(2,2.5))
#'                  )
#'
#' # making 3d plot of joint distribution
#' # require(rgl)
#' rgl::persp3d(res$F1.hat$T1, res$F2.hat$T2, res$F12.hat,
#' xlab = 'T1', ylab = 'T2', zlab = 'F(T1, T2)',
#' main='Joint CDF', col = "lightblue")
#'
#' # add points
#' if (!is.null(res$Pred.Probs)) {
#'     rgl::points3d(res$Pred.Probs$T1, res$Pred.Probs$T2,
#'     res$Pred.Probs$Est.F12, col = "red", add = TRUE)
#' }
#'
#' ## If we want to conduct test for the correlation between two outcomes:
#' \dontrun{
#' if (Sys.info()["sysname"] == "Windows" ){
#'   require(doParallel)
#'   registerDoParallel(5)
#' } else {
#'   require(doMC)
#'   registerDoMC(5)
#' }
#'
#' res = BiIntCensd(SampleDat_case2,
#'                  Corr.Test = TRUE)
#'}
#' @import splines2 CVXR utils grDevices graphics rgl foreach
#' @importFrom stats quantile runif rbinom rnorm pnorm sd
#' @export

BiIntCensd <- function(dat,
                       l = 4,
                       pred.times = NULL,
                       int.k_1 = NULL,
                       int.k_2 = NULL,
                       t1_l = NULL,
                       t1_h = NULL,
                       t2_l = NULL,
                       t2_h = NULL,
                       Corr.Test = FALSE,
                       nBootstrp = 100,
                       seed = 12345,
                       parallel = TRUE) {

  if (length(dat)==2) {
    T1  = dat[[1]]
    T2  = dat[[2]]
  } else if (length(dat)<2){
    stop("No enough data in `dat`! Please check the Help document.")
  } else {
    warning("More than two elements in `dat`, only first two are used!")
  }

  T1 = dat_format(T1, 1)
  T2 = dat_format(T2, 2)
  sample.size = nrow(T1)
  knot = get_knots(T1, T2, int.k_1=int.k_1, int.k_2=int.k_2)
  p_n  = length(knot$knot1)+l
  q_n  = length(knot$knot2)+l


  res = OptLogLik(T1, T2, l = l, knot = knot, p_n = p_n, q_n = q_n,
                  lower1 = t1_l, lower2 = t2_l, upper1 = t1_h, upper2 = t2_h)
  seive.M_ij = res$seive.M_ij
  seive.w_i  = res$seive.w_i
  seive.p_j  = res$seive.p_j

  MC.rho = NULL
  if (Corr.Test) {
    if (parallel) { # parallel computing
      cat("Bootstrapping... \n")
      # progress.bar <- txtProgressBar(0, nBootstrp, style = 3)
      MC.rho = foreach(i = icount(nBootstrp), .combine = c) %dopar% {
        while (TRUE) {
          boot.id = sample(sample.size, replace = T)
          B.res   = suppressWarnings(tryCatch(OptLogLik(T1[boot.id,], T2[boot.id,], l = l, knot = knot, p_n = p_n, q_n = q_n,
                                                        lower1 = t1_l, lower2 = t2_l, upper1 = t1_h, upper2 = t2_h),
                                              error=function(err) {NULL}))
          if (!is.null(B.res)) {
            break
          }
        }
        # setTxtProgressBar(progress.bar, i)
        return(B.res$rho.hat)
      }
    } else { # normal for loop
      cat("Bootstrapping: \n")
      iter = 1
      progress.bar <- txtProgressBar(0, nBootstrp, style = 3)
      while (iter <= nBootstrp) {
        setTxtProgressBar(progress.bar, iter)
        boot.id = sample(sample.size, replace = T)
        B.res   = suppressWarnings(tryCatch(OptLogLik(T1[boot.id,], T2[boot.id,], l = l, knot = knot, p_n = p_n, q_n = q_n,
                                                      lower1 = t1_l, lower2 = t2_l, upper1 = t1_h, upper2 = t2_h),
                                            error=function(err) {NULL}))
        MC.rho  = c(MC.rho, B.res$rho.hat)
        if (!is.null(B.res)) {
          iter    = iter+1
        }
      }
    }

  }

  # estimated F1, F2, and F12.
  T1.seq = sort(c(pred.times[,1], seq(knot$boundary1[1], knot$boundary1[2], length.out = 50)))
  T2.seq = sort(c(pred.times[,2], seq(knot$boundary2[1], knot$boundary2[2], length.out = 50)))
  i.T1   = iSpline(T1.seq, knots=knot$knot1, Boundary.knots=knot$boundary1, degree=l-2, intercept=T)
  i.T2   = iSpline(T2.seq, knots=knot$knot2, Boundary.knots=knot$boundary2, degree=l-2, intercept=T)
  i.T1.col = matrix(rep(i.T1,each=length(T2.seq)),ncol=dim(i.T1)[2], byrow = F)
  i.T2.row = apply(i.T2, 2, rep, length(T1.seq))
  F12.hat  = sapply(seq(nrow(i.T1.col)), function(x) sum(tcrossprod(i.T1.col[x,], i.T2.row[x,]) * seive.M_ij)) # diag(i.T1.col %*% seive.M_ij %*% t(i.T2.row));
  F12.hat = matrix(F12.hat, length(T1.seq), length(T2.seq), byrow = T)
  F1.hat   = i.T1%*%(seive.M_ij%*%matrix(1,q_n-1,1) + seive.w_i)
  F2.hat   = i.T2%*%(t(seive.M_ij)%*%matrix(1,p_n-1,1) + seive.p_j)

  # predicted [F1, F2, F12] for provided times
  if (is.null(pred.times)) {
    Pred.Probs = NULL
  } else {
    # there are better way to do so, but considering that T1/T2 could be identical in multiple tuples [T1, T2], we take a safer way here
    Pred.Probs = NULL
    for (i in seq(nrow(pred.times))) {
      Pred.Probs = rbind(Pred.Probs, c(F1.hat[T1.seq %in% pred.times[i,1]],
                                       F2.hat[T2.seq %in% pred.times[i,2]],
                                       F12.hat[T1.seq %in% pred.times[i,1], T2.seq %in% pred.times[i,2]]
                                       )
                         )
    }
    Pred.Probs = cbind(pred.times, Pred.Probs)
    colnames(Pred.Probs) <- c("T1", "T2", "Est.F1", "Est.F2", "Est.F12")
  }

  cat("Done!\n")
  if (Corr.Test) {
    cat("\nCorrelation test statistic is", res$rho.hat/sd(MC.rho), "corresponding to a p-value:", 1- abs(pnorm(res$rho.hat/sd(MC.rho))-0.5)*2)
  }

  dev.new()
  layout(matrix(c(1,1,1, 2,2,3,3), nrow = 1))
  persp(T1.seq, T2.seq, F12.hat, xlab = 'T1', ylab = 'T2', zlab = 'F(T1, T2)', main='Joint CDF', theta = -35, phi = 20, ticktype = 'detailed')
  plot(T1.seq, F1.hat, col = "black",ylim = c(0,1), type = "l", lwd = 2,  xlab = "T1", ylab = "F(T1)", main = 'Marginal CDF of event 1')
  plot(T2.seq, F2.hat, col = "black",ylim = c(0,1), type = "l",lwd = 2, xlab = "T1", ylab = "F(T2)", main = 'Marginal CDF of event 2')

  return(list(rho.hat    = res$rho.hat,
              MC.rho     = MC.rho,
              seive.M_ij =seive.M_ij,
              seive.w_i  =seive.w_i,
              seive.p_j  =seive.p_j,
              F12.hat    = F12.hat,
              F1.hat     = data.frame(T1 = T1.seq, F1.hat = F1.hat),
              F2.hat     = data.frame(T2 = T2.seq, F2.hat = F2.hat),
              Pred.Probs = data.frame(Pred.Probs)
              )
         )
}




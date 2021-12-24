#' @title Format the input data
#' @param Ti Data from one outcome
#' @param ni Naming the outcome number, i.e., return the column names as tu/tv_<ni>. If not provided, it takes the name from input data.
#' @details See Detail section for \code{BiIntCensd}
#' @return A formatted data with censoring type indicators
#' @export

dat_format <- function(Ti, ni = NULL) {
  if (is.null(ni)) {
    ni = strsplit(colnames(Ti)[1], "_")[[1]][2]
  }
  Ti <- as.data.frame(Ti)
  if (ncol(Ti) == 3) {
    colnames(Ti) <- paste0(c("tu_", "tv_", "Ind_"), ni)
  } else if (ncol(Ti) == 2) {
    ind = NULL
    for (i in seq(nrow(Ti))) {
      ind = c(ind, ifelse(Ti[i,1] == 0, "Left", ifelse(is.infinite(Ti[i,2]), "Right", "Interval")))
    }
    Ti = cbind(Ti, ind)
    colnames(Ti) <- paste0(c("tu_", "tv_", "Ind_"), ni)
  } else {
    stop("The T1 and T2 in `dat` is not correctly formatted!")
  }
  return(Ti)
}

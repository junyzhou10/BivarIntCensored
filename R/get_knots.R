#' @title Determine the position of knots
#' @param T1/T2 Observation intervals
#' @param int.k_1/int.k_2 Prespicified number of interval knots
#' @details If int.k_1/2 not given, then the # of interior knots = \eqn{n^(1/3)}; Otherwise, equals to int.k_1/2
#' @return A list of two elements, each contains the position of interval knots
#' @noRd

get_knots<-function(T1, T2, int.k_1 = NULL, int.k_2 = NULL){
  if ( is.null(int.k_1) ) {
    int.k_1 = round(dim(T1)[1]^(1/3))
  }

  if ( is.null(int.k_2) ) {
    int.k_2 = round(dim(T2)[1]^(1/3))
  }

  knot_pre_1 = c(T1$tu_1, T1$tv_1)
  knot_pre_1 = knot_pre_1[is.finite(knot_pre_1)]
  knot_pre_1 = knot_pre_1[knot_pre_1>0]
  knot_pre_2 = c(T2$tu_2, T2$tv_2)
  knot_pre_2 = knot_pre_2[is.finite(knot_pre_2)]
  knot_pre_2 = knot_pre_2[knot_pre_2>0]

  pos1 = quantile(knot_pre_1, seq(0,1,length.out = int.k_1+2), names = F)
  pos2 = quantile(knot_pre_2, seq(0,1,length.out = int.k_2+2), names = F)
  return(list(knot1=pos1[2:(length(pos1)-1)],
              knot2=pos2[2:(length(pos2)-1)],
              boundary1 = c(max(0, pos1[1]-0.5), pos1[length(pos1)]+0.5),
              boundary2 = c(max(0, pos2[1]-0.5), pos2[length(pos2)]+0.5)
  ))
}

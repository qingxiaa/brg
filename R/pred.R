#' @title pred
#' @description FUNCTION_DESCRIPTION
#' @param alpha_EB PARAM_DESCRIPTION
#' @param beta_EB PARAM_DESCRIPTION
#' @param bridge_ind PARAM_DESCRIPTION
#' @param Z1t1 PARAM_DESCRIPTION
#' @param Z2t1_M PARAM_DESCRIPTION
#' @param J PARAM_DESCRIPTION
#' @param M PARAM_DESCRIPTION
#' @param N PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @noRd
pred <- function(alpha_EB, beta_EB, bridge_ind, Z1t1, Z2t1_M, J, M, N){
  Z2t1_P <- matrix(NA, N, J)
  for ( j in 1:J){
    Z2t1_P[bridge_ind, j] <- Z2t1_M[, j]
    
    for (q in 1:(N-M)){
      Z2t1_P[setdiff(1:N, bridge_ind)[q], j] <- (alpha_EB[j]  + Z1t1[setdiff(1:N, bridge_ind)[q], j] * beta_EB[j])
    }
  }
  return(Z2t1_P)
}
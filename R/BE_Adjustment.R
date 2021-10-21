#' @title BE_Adjustment
#' @description FUNCTION_DESCRIPTION
#' @param phi_EB PARAM_DESCRIPTION
#' @param delta_EB PARAM_DESCRIPTION
#' @param Z PARAM_DESCRIPTION
#' @param mu_hat PARAM_DESCRIPTION
#' @param sigma_hat PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @noRd
BE_Adjustment <- function(phi_EB, delta_EB, Z, mu_hat, sigma_hat){
  adj_y <- sigma_hat/delta_EB * (Z - phi_EB) + mu_hat
  return(adj_y)
}
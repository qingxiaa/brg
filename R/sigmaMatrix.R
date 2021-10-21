
#' @title sigmaMatrix
#' @description FUNCTION_DESCRIPTION
#' @param delta1j PARAM_DESCRIPTION
#' @param delta2j PARAM_DESCRIPTION
#' @param rho1 PARAM_DESCRIPTION
#' @param rho2 PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @noRd
sigmaMatrix <- function(delta1j, delta2j, rho1, rho2) {
  rho3 <- rho1 * rho2
  Sigmaj <- matrix(c(delta1j^2, rho1 * delta1j * delta2j, rho1 * rho2 * delta1j * delta2j, 
                     rho1 * delta1j * delta2j, delta2j^2, rho2 * delta2j^2, 
                     rho1 * rho2 * delta1j * delta2j, rho2 * delta2j^2, delta2j^2), 3, 3)
  return(Sigmaj)
}

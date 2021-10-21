#' @title sigmaUpdate
#' @description FUNCTION_DESCRIPTION
#' @param lambda0 PARAM_DESCRIPTION
#' @param theta0 PARAM_DESCRIPTION
#' @param betaVec PARAM_DESCRIPTION
#' @param X PARAM_DESCRIPTION
#' @param y PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[MCMCpack]{InvGamma}}
#' @noRd
#' @importFrom MCMCpack rinvgamma
sigmaUpdate <- function(lambda0, theta0, betaVec, X, y){
  n <- length(y)
  ssE <- SSE(betaVec, X, y)
  sigmaUpd <- MCMCpack::rinvgamma(1, lambda0  + n/2, theta0  + ssE/2)
}
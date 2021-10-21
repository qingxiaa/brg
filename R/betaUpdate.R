#' @title betaUpdate
#' @description FUNCTION_DESCRIPTION
#' @param sigma PARAM_DESCRIPTION
#' @param beta0 PARAM_DESCRIPTION
#' @param X PARAM_DESCRIPTION
#' @param y PARAM_DESCRIPTION
#' @param precS PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[MASS]{mvrnorm}}
#' @noRd
#' @importFrom MASS mvrnorm
betaUpdate <- function(sigma, beta0, X, y, precS){ #sigma_square
  sigma_beta_prec <- precS/sigma # it is squared sigma
  meanVec <- solve(sigma_beta_prec + t(X) %*% X/sigma) %*% (sigma_beta_prec %*% beta0 + t(X) %*% y/sigma)
  VarMat <- solve(sigma_beta_prec + t(X) %*% X/sigma)
  betaVec <- MASS::mvrnorm(n = 1, meanVec, VarMat)
  return(betaVec)
}
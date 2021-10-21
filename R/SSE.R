#' @title SSE
#' @description FUNCTION_DESCRIPTION
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
#' @noRd
SSE <- function(betaVec, X, y){
  out <- t(y - X %*% betaVec) %*% (y - X %*% betaVec)
  return(out)
}
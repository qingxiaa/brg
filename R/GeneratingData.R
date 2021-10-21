# Generate data function
#' @title GeneratingData
#' @description FUNCTION_DESCRIPTION
#' @param N PARAM_DESCRIPTION
#' @param J PARAM_DESCRIPTION
#' @param mu01j PARAM_DESCRIPTION
#' @param mu02j PARAM_DESCRIPTION
#' @param mu PARAM_DESCRIPTION
#' @param Sigma PARAM_DESCRIPTION
#' @param sigma01j PARAM_DESCRIPTION
#' @param sigma02j PARAM_DESCRIPTION
#' @param phi1j PARAM_DESCRIPTION
#' @param phi2j PARAM_DESCRIPTION
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
GeneratingData=function(N, J, mu01j, mu02j, mu, Sigma, sigma01j, sigma02j, phi1j, phi2j){
  data <- matrix(NA, N, J)
  Z1_t1= matrix(NA, N, J)
  Z2_t1= matrix(NA, N, J)
  Z2_t2= matrix(NA, N, J)
  Y1_t1= matrix(NA, N, J)
  Y2_t1= matrix(NA, N, J)
  Y2_t2= matrix(NA, N, J)
  for (j in 1:J)
  {
    data <- MASS::mvrnorm(N, mu[, j], Sigma[[j]])
    Z1_t1[, j] <- data[, 1]
    Z2_t1[, j] <- data[, 2]
    Z2_t2[, j] <- data[, 3]
    Y1_t1[, j] <- (Z1_t1[, j] - phi1j[j]) * sigma01j[j] + phi1j[j] + mu01j[j]
    Y2_t1[, j] <- (Z2_t1[, j] - phi2j[j]) * sigma01j[j] + phi2j[j] + mu01j[j]
    Y2_t2[, j] <- (Z2_t2[, j] - phi2j[j]) * sigma02j[j] + phi2j[j] + mu02j[j]
  }
  mylist <- list(Y1_t1, Y2_t1, Y2_t2)
  return(mylist)
}

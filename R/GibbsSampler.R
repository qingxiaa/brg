#' @title GibbsSampler
#' @description Gibbs Sampler for posterior sampling in BRIDGE
#' @param y PARAM_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param beta0 PARAM_DESCRIPTION
#' @param precS PARAM_DESCRIPTION
#' @param lambda0 PARAM_DESCRIPTION
#' @param theta0 PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @noRd
GibbsSampler <- function(y, x, beta0, precS, lambda0, theta0){
  n <-length(y)
  X <- cbind(rep(1, n), x)
  
  y=y
  n_params <- 2 
  n_iterations <- 500 
  burn_in=150
  
  
  betaVec <- c(0, 0) # starting value
  sigmaUpd <- 0.1 #starting value
  
  beta_out <- matrix(data = NA, nrow = n_iterations, ncol = n_params)
  sigma_out <- matrix(data <- NA, nrow = n_iterations, ncol = 1)
  for (i in 1:n_iterations){
    betaVec <- betaUpdate(sigmaUpd, beta0, X, y, precS)
    sigmaUpd <- sigmaUpdate(lambda0, theta0, betaVec, X, y)
    beta_out[i, ] <- betaVec
    sigma_out[i, ] <- sigmaUpd
  }
  out <- list(beta_out[(burn_in + 1):n_iterations, ], sigma_out[(burn_in + 1):n_iterations, ])
  return(out) # out is a list if length <- J, in j list, a nested list with c(alpha, beta) and sigma^2
}


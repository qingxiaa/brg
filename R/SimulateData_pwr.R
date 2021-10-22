#'@title This function is to simulate genes expression data for evaluating the performance of BRIDGE
#'@description This function is to simulate genes expression data for evaluating the performance of BRIDGE
#'@param J The number of genes
#'@param N The number of subjects at each time point
#'@param G The number of differential genes
#'@param M The number of bridging samples
#'@param TimeEff0 Time effects across two time points for differential genes
#'@param mu01 mean expression at time point 1
#'@param sigma01 variance of gene expression at time point 1
#'@param sigma02 variance of gene expression at time point 2
#'@param delta2 multiplicative batch effects of batch 2
#'@param rho1 correlation coefficient between different batch
#'@param rho2 correlation coefficient between different time point 
#'@param phi2 additive batch effects of batch 2
#'@return Function outputs a list including the following:
#' \describe{
#'     \item{\code{Data}}{A simulated genomic measure matrix (observation x genes) with batch, time and subjectID information at first 3 columns }
#'     \item{\code{bridgingSampleInd}}{The index of subjects that are bridging samples}
#'     \item{\code{Geneind}}{The index of genes that are truely differential}
#'     \item{\code{Obs_te}}{The vector of observed time effect size of truly differential genes}}
#' @import MASS 
#' @noRd
SimulateData_pwr <- function(J, N, G, M, TimeEff0, mu01, sigma01, sigma02, delta2, rho1, rho2, phi2){
  phi1j <- phi2j <- NULL
  delta1 <- 1 # batch 1 is reference batch
  phi1 <- 0 # batch 1 is reference batch
  
  # sample parameters for each gene
  a <- 0.1
  b <- 0.1
  delta1j <- runif(J, delta1 - a, delta1 + a)
  delta2j <- runif(J, delta2 - a, delta2 + a)
  sigma01j <- runif(J, sigma01 - b, sigma01 + b)
  sigma02j <- runif(J, sigma02 - b, sigma02 + b)
  phi1j <- rnorm(J, phi1, a)
  phi2j <- rnorm(J, phi2, a)
  mu01j <- rnorm(J, mu01, b)
  Geneind <- sample(J, G)
  Ifun <- (1:J) %in% Geneind
  
  te0 <- rnorm(J, TimeEff0, 0.1)
  # Get mu02j, mu and Sigma to call GeneratingData()
  mu02j <- mu01j  + Ifun * te0
  
  Obs_te <- te0[Geneind]
  
  mu <- matrix(NA, 3, J)
  for (j in 1:J){
    mu[, j] <- c(phi1j[j], phi2j[j], phi2j[j])
  }
  
  Sigma  <- list()
  for (j in 1:J)
  {
    Sigma[[j]] <- sigmaMatrix(delta1j[j], delta2j[j], rho1, rho2)
  }
  Data <- GeneratingData(N, J, mu01j, mu02j, mu, Sigma, sigma01j, sigma02j, phi1j, phi2j)
  
  Y1_t1 <- Data[[1]] 
  Y2_t1 <- Data[[2]]
  Y2_t2 <- Data[[3]]
  
  sample_ind <- sample(c(1:N), M)
  Y2_t1_M <- Y2_t1[sample_ind, ]
  
  batch <- c(rep("batch1",N),rep("batch2",M+N))
  time <- c(rep("time1",N+M),rep("time2",N))
  subjectId <- c(1:N,sample_ind,1:N)
  Data <- rbind(Y1_t1, Y2_t1_M, Y2_t2)
  df <- data.frame(batch, time, subjectId, Data)
  out <- list(Data = df, 
              bridgingSampleInd = sample_ind, 
              Geneind = Geneind,
              Obs_te = Obs_te)
  return(out)
}

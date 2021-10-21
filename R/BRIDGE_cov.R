#' @title Batch effect Reduction of mIcroarray data with Dependent samples usinG an Empirical Bayes approach (BRIDGE) - Covariates. 
#' @description  To circumvent the limitations of existing batch-effect correction methods studies involving the longitudinal collection of high-dimensional 'omic' data, we propose Batch effect coRrectIon of microarray data with Dependent samples usinG an Empirical Bayes approach (BRIDGE). BRIDGE accounts for within-subject dependency expected in longitudinal studies with high-dimensional 'omic' data and involves the estimation and correction of additive batch and multiplicative batch effects when batch effects are confounded with time. After correcting for batch effects with BRIDGE, adjusted data can be used for downstream statistical analyses as if all samples are measured in a single batch. BRIDGE is applicable to different 'omics' platforms, such as microarray gene expression, DNA methylation or neuroimaging studies of neurodevelopment, as long as the transformed data are high-dimension and approximately normal distribution.
#'@param infile A genomic measure matrix (genes x observations) - for example, gene expression matrix.  \code{NA}s are not allowed
#'@param SubjectID A vector of the subject ID of the 'infile' data matrix. The samples are collected from the same subject having a unique ID.   
#'@param batch A factor vector indicates each observation of the 'infile' data matrix are obtained from which batches. The labels of batch are \code{batch1} and \code{batch2}. 
#'@param time A factor vector indicates each observation of the 'infile' data matrix are collected at which time points. The labels of time points are \code{time1} and \code{time2}.
#'@param mdl A A design matrix including covariates information.
#'@param data a data frame containing all information except expression data (observation x variables) including batch, time, SubjectID and other covariates.
#'@param optimization (Optical) Default value is \code{FALSE}. If \code{TRUE}, BRIDGE automatically choose one method of \code{imputation} to adjusts data in order to achieve higher power (maybe) at the expense of inflated type I error. 
#'@param imputation (Optical) Default value is \code{FALSE}. When \code{optimization = FALSE}, \code{imputation = FALSE} adjust data using empirical Bayes estimations of additive batch effects and multiplicative batch effects and after correcting the batch effects(Method 1), the adjusted data can be used for downstream analysis as if all samples are measured in batch 1; and \code{imputation = TRUE}, BRIDGE adjusts data based on linear regression imputation using empirical Bayes estimate of linear coefficients for the date of the non-bridging samples(Method 2). Method 2 works better in controlling type I error. After correcting the batch effects, the adjusted data can be used for downstream analysis as if all samples are measured in batch 2.
#'@return Function outputs a list including the following:
#' \describe{
#'     \item{\code{factor.df}}{A data frame contains batch, time and subject id information}
#'     \item{\code{AdjustedData}}{A genomic measure matrix (genes x observations) with batch, time and subjectID information at first 3 rows. If method 2 is used, then return data matrix adds new columns containing imputed data for non-bridging samples.}
#'     \item{\code{phi2hat}}{A vector of J genes/features containing empirical Bayes estimate of additive batch effects og batch 2 for each feature, assuming phi1 = 0}
#'     \item{\code{delta2hat}}{A vector of J genes/features containing empirical Bayes estimate of multiplicative batch effects og batch 2 for each feature, assuming delta1 = 1}
#'     \item{\code{alphahat}}{A vector of J genes/features containing empirical Bayes estimate of linear coefficient}
#'     \item{\code{betahat}}{A vector of J genes/features containing empirical Bayes estimate of linear coefficient}
#'     }
#'@import MASS MCMCpack stats
#'@export
#'@rdname BRIDGE_cov


BRIDGE_cov <- function(infile, SubjectID, batch, time, mdl=NULL, data, imputation = FALSE, optimization = FALSE){
  #############Step 0: Prepare the Data###################
  dataMat <- as.matrix(t(infile)) # observation x genes matrix
  # check for missing data 
  if (sum(is.na(dataMat)) > 0) {
    message <- paste0('Missing data found before running BRIDGE')
    stop(message)
  }
  # number of observation 
  ObsSize <- nrow(dataMat) 
  Obs <- c(1:ObsSize)
  # extract the levels of batch and time, assume the levels of batch or time are in chronological sequence 1,2...
  batch1 <- names(sort(table(batch), decreasing <- FALSE)[1])
  batch2 <- names(sort(table(batch), decreasing <- FALSE)[2])
  time1 <- names(sort(table(time), decreasing <- TRUE)[1])
  time2 <- names(sort(table(time), decreasing <- TRUE)[2])
  
  # Extract observation Id, subject Id of bridging samples and size of bridging samples
  bridge_Obs <- which(batch == batch2 & time == time1)
  bridge_id <- SubjectID[bridge_Obs]
  
  # Extract other observation Id's
  Y1t1_Obs <- which(batch == batch1 & time == time1)
  Y1t1_id <- SubjectID[Y1t1_Obs]
  Y2t2_Obs <- which(batch == batch2 & time == time2)
  Y2t2_id <- SubjectID[Y2t2_Obs]
  
  # number of subjects
  N1 <- length(Y1t1_Obs)
  N2 <- length(Y2t2_Obs)
  
  # number of genes
  J <- ncol(dataMat)
  
  # confirm bridge_id sample in batch 1
  bridge_id <- bridge_id[bridge_id %in% Y1t1_id]
  M <- length(bridge_id)
  # prepare data for BRIDGE
  Y1_t1 <- dataMat[Y1t1_Obs,]
  bridge_ind <- match(bridge_id, Y1t1_id) # the index(integer) of bridging samples in Y1_t1 
  Y1_t1_M <- Y1_t1[bridge_ind,]
  
  Y2_t1_M <- dataMat[bridge_Obs,][SubjectID[bridge_Obs] %in% bridge_id, ]
  Y2_t2 <- dataMat[Y2t2_Obs,]
  
  
  
  
  
  if (N1 == N2){
    cat("[BRIDGE] found", ObsSize,"total observations", N1, "subjects", J, "features\n", M, "bridging samples")
  }else{
    cat("[BRIDGE] found", ObsSize,"total observations", N1, 
        "subjects in time 1 and ",N2, "subjects in time 2,", J, "features\n", M, "bridging samples")
  }
  ## Method for batch effects correction###
  
  ############ Step 1: Data normalized ########################
  
  # Get all estimates of parameters to normalize Y to Z
  batchmod <- model.matrix(~batch)  
  timemod <- model.matrix(~-1+time)  
  
  
  #ref <- which(levels(as.factor(batch))==batch1) # find the reference
  
  design <- cbind(batchmod,timemod,mdl)
  
  check <- apply(design, 2, function(x) all(x == 1))
  #throw away the intercept
  design <- as.matrix(design[,!check])
  
  beta_hat <- solve(crossprod(design), tcrossprod(t(design), t(dataMat)))
  grand.design <- design
  grand.design[,1] <- 0
  grand.means <- grand.design %*% beta_hat
  
  pred.means <- design %*% beta_hat
  
  bias <- dataMat - pred.means
  
  sigma01j_h <- apply(bias[Y1t1_Obs,],2,sd)
  delta2j_h <- apply(bias[bridge_Obs,],2,sd)/sigma01j_h
  sigma02j_h <- apply(bias[Y2t2_Obs,],2,sd)/delta2j_h
  
  sigma_h <- t(cbind(sigma01j_h %*% t(rep(1, N1+M)),sigma02j_h %*% t(rep(1, N2))))
  
  
  
  # Normalize data
  cat("[BRIDGE] standardizing data across features\n")
  
  Z.data <- bias/sigma_h + (pred.means-grand.means)
  
  Z1t1 <- Z.data[Y1t1_Obs,]
  Z1t1_M <- Z1t1[bridge_ind,]
  Z2t1_M <- Z.data[bridge_Obs,][SubjectID[bridge_Obs] %in% bridge_id, ]
  Z2t2 <- Z.data[Y2t2_Obs,]
  
  
  
  
  ########## Step 2: Empirical Bayes Estimate of Hyperprior and LS parameters########
  md <- vector(mode <- 'list', length <- J) # store models
  
  sigma2_hat <- NULL
  a0_hat <- NULL 
  b0_hat <- NULL
  Sa_hat <- NULL
  Sb_hat <- NULL
  Sab_hat <- NULL
  SSE <- NULL
  
  # Fit series of ordinary linear regression
  cat(paste0('[BRIDGE] fitting lm model for features \n'))
  
  
  for ( j in 1:J){  
    md[[j]] <- lm(Z2t1_M[, j] ~ Z1t1_M[, j])
    sigma2_hat[j] <- sum(md[[j]]$residuals^2)/(M-1)
    a0_hat[j] <- md[[j]]$coefficients[1]
    b0_hat[j] <- md[[j]]$coefficients[2]
    Sa_hat[j] <- vcov(md[[j]])[1, 1]/sigma2_hat[j]
    Sb_hat[j] <- vcov(md[[j]])[2, 2]/sigma2_hat[j]
    Sab_hat[j] <- vcov(md[[j]])[1, 2]/sigma2_hat[j]
  }
  
  # beta normal prior
  a0_bar <- mean(a0_hat)
  b0_bar <- mean(b0_hat)
  beta0 <- c(a0_bar, b0_bar)
  # precison matrix
  Sa_bar <- mean(Sa_hat)
  Sb_bar <- mean(Sb_hat)
  Sab_bar <- mean(Sab_hat)
  precS <- solve(matrix(c(Sa_bar, Sab_bar, Sab_bar, Sb_bar), 2, 2))
  # tau inverse gamma prior 
  V_bar <-  mean(sigma2_hat)
  S2_bar <- var(sigma2_hat)
  lambda0 <- V_bar^2/S2_bar + 2
  theta0 <- V_bar^3/S2_bar + V_bar
  
  # 2) Gibbs sampling for EB LS parameter estimates
  Gibbs <- list()
  Gibbs <- lapply(1:J,function(x) GibbsSampler(Z2t1_M[,x],Z1t1_M[,x], beta0, precS, lambda0, theta0))
  sigma_EB <- sapply(1:J,function(x) mean(sqrt(Gibbs[[x]][[2]])))
  alpha_EB <- sapply(1:J,function(x) apply(Gibbs[[x]][[1]],2,mean))[1,]
  beta_EB <- sapply(1:J,function(x) apply(Gibbs[[x]][[1]],2,mean))[2,]
  
  
  # Empirical Bayes Estimates for method 1 
  phi2_EB <- alpha_EB
  delta2_EB <- sqrt(beta_EB^2 + sigma_EB^2)
  
  ########################## Step 3: ADJUST DATA######################
  
  # Adjust data using Method 1
  if ( (imputation == FALSE & optimization == FALSE) | (optimization == TRUE & mean(delta2_EB) > 1)){
    cat("[BRIDGE] method 1 is used. \n")
    
    deltaMat <- rbind(matrix(1,N1,J),t(delta2_EB %*% t(rep(1,N2+M))))
    phiMat <- rbind(matrix(0,N1,J),t(phi2_EB %*% t(rep(1,N2+M))))
    adj_DataMat <-sigma_h/deltaMat * (Z.data - phiMat) + grand.means
    factor.df <- data
    AdjustedData <- t(adj_DataMat)
  } 
  else if( (imputation == TRUE & optimization == FALSE) | (optimization == TRUE & mean(delta2_EB) <= 1)){     # Adjust data using method 2 based on linear regression imputation   
    cat("[BRIDGE] method 2 is used. \n")
    
    Z2t1_P <- t(sapply(1:N1, function(x) t(alpha_EB  + Z1t1[x,] * beta_EB)))
    Z2t1_P[bridge_ind,] <- as.matrix(Z2t1_M)
    
    deltaMat <- t(delta2_EB %*% t(rep(1,N1)))
    phiMat <- t(phi2_EB %*% t(rep(1,N1)))
    Y2t1_P <- sigma_h[Y1t1_Obs,] * (Z2t1_P - phiMat) + phiMat + grand.means[Y1t1_Obs,]
    
    adj_DataMat <- rbind(Y1_t1, Y2t1_P, Y2_t2)
    labs <- rep("Imputed",N1)
    labs[bridge_ind] = batch2
    imputed.label <- labs
    new.batch <- c(rep(batch1, N1),imputed.label,rep(batch2,N2))
    factor.df <- rbind(data[Y1t1_Obs,],data[Y1t1_Obs,],data[Y2t2_Obs,])
    wh.is.batch <- which(colnames(factor.df)=="batch")
    factor.df [,wh.is.batch] <- new.batch
    AdjustedData <- t(adj_DataMat)
  }
  
  
  
  ################ Finished ##################
  return(list(factor.df = factor.df,
              AdjustedData = AdjustedData,
              phi2hat = phi2_EB,
              delta2hat = delta2_EB,
              alphahat = alpha_EB,
              betahat = beta_EB))
}

library(faraway)

"%+%" <- function(e1, e2) {
  paste(e1, e2, sep = "")
}

create_logger <- function(log_file) {
  force(log_file)
  
  function(...) {
    txt <- paste(...)
    cat(txt, "\n")                                   # to console
    cat(txt, "\n", file = log_file, append = TRUE)   # to file
  }
}

# Helper: build Sigma from (a1,a2,a3)
make_Sigma <- function(a1, a2, a3) {
  tau1 <- exp(a1)
  tau2 <- exp(a2)
  rho <- tanh(a3)
  matrix(c(
    tau1^2, rho * tau1 * tau2,
    rho * tau1 * tau2, tau2^2
  ), 2, 2, byrow = TRUE)
}

# Helper: stable chol-based solver with jitter
# Adds jitter * I to M and attempts Cholesky; if it still fails, sets fail = TRUE.
logdet_and_solve <- function(M, jitter = 1e-10) {
  Mj <- M + diag(jitter, 2)
  dec <- tryCatch(chol(Mj), error = function(e) NULL)
  if (is.null(dec)) {
    return(list(fail = TRUE))
  }
  logdet <- 2 * sum(log(diag(dec)))
  invM <- chol2inv(dec)
  list(fail = FALSE, logdet = logdet, inv = invM, chol = dec)
}

# Compute the observed logit sensitivities and logit specificities
comp.y <- function(Data){
  k <- nrow(Data)
  # A continuity correction of 0.5 will be added to a cell which appears to contain 0, since logit(p) won't be defined otherwise
  
  y <- matrix(c(logit(madad(Data, correction.control = "single")$sens$sens), 
                logit(madad(Data, correction.control = "single")$spec$spec)), nrow=k, ncol=2)
  return(y)
}

# Compute the within-study variances of logit Se and logit Sp using the delta-method                                                                           #
comp.Psi <- function(Data){
  # "Data" is the kby4 data frame containing the number of "TP", "TN", "FP" and "FN".
  k <- nrow(Data); p <- 2
  correction <- 0.5
  Y <- comp.y(Data=Data)
  
  Psi <- list()
  for(i in 1:k){
    if(Data$TP[i]==0|Data$FP[i]==0|Data$FN[i]==0|Data$TN[i]==0){
      Data$TP[i] = Data$TP[i] + correction
      Data$FP[i] = Data$FP[i] + correction
      Data$FN[i] = Data$FN[i] + correction
      Data$TN[i] = Data$TN[i] + correction
    }
    n1 <- Data$TP+Data$FN; n2 <- Data$FP+Data$TN
    Psi[[i]] <- matrix(c((1/(n1[i]*(ilogit(Y[i,1])*(1-ilogit(Y[i,1]))))), 0, 0, (1/(n2[i]*(ilogit(Y[i,2])*(1-ilogit(Y[i,2])))))), 2,2)
  }
  
  return(Psi)
}
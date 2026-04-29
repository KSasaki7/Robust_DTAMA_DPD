source("./src/functions/utils.R")
library(mada)
library(numDeriv)

# ---- Robust DPD-based method under the Bivariate normal-normal model (BNN-DPD) ----
# Data inputs are assumed to be in the form of 2x2 tables (TP, FN, FP, TN) for each study.
# alpha > 0 : DPD tuning parameter (alpha = 0 -> MLE fallback)
fit_bnnDPD <- function(dat, alpha, verbose=FALSE){
  df <- make_logit_data(TP=dat$TP, FN=dat$FN, FP=dat$FP, TN=dat$TN, 
                        correction.control = "single")
  
  # --- initial values from standard BNN (alpha = 0) ---
  res_bnn = fit_bnn(Data = dat, method = "ml")
  
  check_mu  <- any(!is.finite(res_bnn$SeSp.logit))
  check_Psi <- any(!is.finite(diag(res_bnn$Psi)))
  if (check_mu) {
    warning("fit_bnn returned non-finite SeSp.logit; default initial values are used.")
  }
  if (check_Psi) {
    warning("fit_bnn returned non-finite diagonal elements of Psi; default initial values are used.")
  }
  
  if (check_mu || check_Psi){
    init.par_trans <- NULL
  } else {
    MLE.par_nat <- c(
      res_bnn$SeSp.logit,                                    # mu1, mu2
      sqrt(res_bnn$Psi[1, 1]),                               # tau1
      sqrt(res_bnn$Psi[2, 2]),                               # tau2
      res_bnn$Psi[1, 2] / sqrt(res_bnn$Psi[1, 1] * res_bnn$Psi[2, 2]) # rho
    )
    names(MLE.par_nat) <- c("mu1", "mu2", "tau1", "tau2", "rho")
    
    tau1_safe <- max(MLE.par_nat["tau1"], 1e-4)
    tau2_safe <- max(MLE.par_nat["tau2"], 1e-4)
    rho_clipped <- min(max(MLE.par_nat["rho"], -0.999), 0.999)
    
    init.par_trans <- c(
      mu1 = MLE.par_nat["mu1"],
      mu2 = MLE.par_nat["mu2"],
      a1  = log(tau1_safe),
      a2  = log(tau2_safe),
      a3  = atanh(rho_clipped)
    )
  }
  
  fit_DPD = fit_bnnDPD_core(df$yi1, df$yi2, df$si1, df$si2, 
                               init.par = init.par_trans,
                               alpha=alpha,
                               verbose = verbose)
  return (fit_DPD)
}

fit_bnnDPD_core <- function(yi1, yi2, si1, si2, alpha,
                           init.par = NULL, 
                           compute_variance = TRUE,
                           seed = NULL,
                           control = list(maxit = 1000, reltol = 1e-4),
                           jitter = 1e-10, verbose = TRUE) {
  stopifnot(
    length(yi1) == length(yi2),
    length(yi1) == length(si1),
    length(yi1) == length(si2),
    alpha >= 0
  )
  control_default <- list(maxit = 1000, reltol = 1e-4)
  control <- modifyList(control_default, control)
  
  n <- length(yi1)
  Y <- cbind(yi1, yi2)
  # Within-study variance matrices are diagonal
  Si_list <- lapply(seq_len(n), function(i) diag(c(si1[i]^2, si2[i]^2), nrow = 2))
  
  # Unconstrained parameterization to ensure Sigma is PD:
  # tau1 = exp(a1), tau2 = exp(a2), rho = tanh(a3), so
  # Sigma = [[tau1^2, rho*tau1*tau2],[rho*tau1*tau2, tau2^2]]
  # mu = (m1, m2)
  if (is.null(init.par)) {
    # crude starting values:
    w1 <- 1 / (si1^2 + median(si1^2))
    w2 <- 1 / (si2^2 + median(si2^2))
    m1_0 <- sum(w1 * yi1) / sum(w1)
    m2_0 <- sum(w2 * yi2) / sum(w2)
    a1_0 <- log(0.2) # tau1 ~ 0.2
    a2_0 <- log(0.2) # tau2 ~ 0.2
    a3_0 <- atanh(0) # rho ~ 0
    init.par <- c(m1_0, m2_0, a1_0, a2_0, a3_0)
  } else if (length(init.par) != 5L) {
    stop("init.par must be a numeric vector of length 5: c(mu1, mu2, a1, a2, a3).")
  }
  
  # D_alpha objective (we maximize it; optim does minimization, so return -D)
  D_alpha_obj <- function(m1, m2){
    function(par) {
      theta <- c(m1, m2, par)
      D_i_vec <- vapply(
        seq_len(n),
        function(i) dpd_contrib_i(theta, i, Y, Si_list, alpha, jitter),
        numeric(1)
      )
      val = -sum(D_i_vec)

      val
    }
  }
  
  update_mu_fixed_point = function(m1,m2,a1,a2,a3,jitter){
    mu <- c(m1, m2)
    Sigma <- make_Sigma(a1, a2, a3)
    c2pi <- log(2 * pi)
    
    gW_sum = matrix(0,nrow = 2, ncol=2)
    gWy_sum = matrix(0, nrow=2, ncol=1)
    
    for (i in 1:n) {
      Vi <- Si_list[[i]] + Sigma
      ld <- logdet_and_solve(Vi, jitter = jitter)
      if (ld$fail) {
        return(rep(NA_real_, 2))
      }
      
      diff <- (Y[i, ] - mu)
      quad <- as.numeric(t(diff) %*% ld$inv %*% diff) # Mahalanobis^2
      logQi <- -c2pi - 0.5 * ld$logdet # log Q_i
      

      # DPD pieces using log-space for stability
      # g_i = exp( alpha*logQi - 0.5*alpha*quad )
      log_gi <- alpha * logQi - 0.5 * alpha * quad
      gi <- exp(log_gi)
      
      giWi = gi * ld$inv
      giWiyi = gi * ld$inv %*% Y[i, ]
      
      gW_sum <- gW_sum + giWi
      gWy_sum <- gWy_sum + giWiyi
      
    }
    ld_gW <- logdet_and_solve(gW_sum, jitter = jitter)
    if (ld_gW$fail) {
      return(rep(NA_real_, 2))
    }
    mu_new = ld_gW$inv %*%  gWy_sum
    return(mu_new)
  }
  
  m1_0 =init.par[1]
  m2_0 =init.par[2]
  a1_0 = init.par[3]
  a2_0 = init.par[4]
  a3_0 = init.par[5]
  bounds = list(a1 = c(log(1e-4), log(10)),
                a2 = c(log(1e-4), log(10)),
                a3 = c(atanh(-0.999), atanh(0.999)))
  lower <- c(bounds$a1[1], bounds$a2[1], bounds$a3[1])
  upper <- c(bounds$a1[2], bounds$a2[2], bounds$a3[2])
  
  opt <- NULL
  success <- TRUE
  
  for (num in 1:control$maxit){
    mu = update_mu_fixed_point(m1_0,m2_0,a1_0,a2_0,a3_0,jitter)
    if (any(!is.finite(mu))) {
      warning("update_mu_fixed_point failed; returning NA")
      success <- FALSE
      break
    }
    m1 = mu[1]
    m2 = mu[2]
    
    opt <- optim(c(a1_0,a2_0,a3_0), D_alpha_obj(m1, m2),
                 method = "L-BFGS-B",
                 lower = lower, upper = upper
    )
    a1 <- opt$par[1]
    a2 <- opt$par[2]
    a3 <- opt$par[3]
    
    diff = c(abs(m1 - m1_0), abs(m2 - m2_0),
             abs(a1 - a1_0), abs(a2 - a2_0), abs(a3 - a3_0))
    rmse <- max(diff)
    
    if (rmse < control$reltol) {
      break
    }
    # update
    m1_0 <- m1
    m2_0 <- m2
    a1_0 <- a1
    a2_0 <- a2
    a3_0 <- a3
    
  }
  
  if (!success || is.null(opt)) {
    par_trans_hat <- c(mu1 = NA_real_, mu2 = NA_real_,
                 a1 = NA_real_, a2 = NA_real_, a3 = NA_real_)
    par <- c(mu1 = NA_real_, mu2 = NA_real_,
             tau1 = NA_real_, tau2 = NA_real_, rho = NA_real_)
    
    return(list(
      est_summary = NA,
      SeSp = c(sensitivity = NA_real_, specificity = NA_real_),
      par = par,
      Sigma = matrix(NA_real_, 2, 2),
      cov = NA,
      alpha = alpha,
      H_score = NA_real_,
      weights_mu = NA,
      par_trans = par_trans_hat,
      opt.obj = NULL,
      value = NA_real_,
      convergence = NA_integer_,
      message = "mu-update failed or optimization was not performed.",
      call = match.call()
    ))
  }
  
  par_trans_hat <- c(m1, m2, a1, a2, a3)
  names(par_trans_hat) <- c("mu1", "mu2", "a1", "a2", "a3")

  if (verbose) {
    cat("Convergence:", opt$convergence, " (0 is good)\n")
    cat("Message:", opt$message, "\n")
  }
  
  # Back-transform estimates
  tau1 <- exp(a1)
  tau2 <- exp(a2)
  rho <- tanh(a3)
  Sigma_hat <- make_Sigma(a1, a2, a3)
  
  # Back-transform pooled sensitivity/specificity
  pooled_sens <- plogis(m1)
  pooled_spec <- plogis(m2)
  par <- setNames(
    c(m1, m2, tau1, tau2, rho),
    c("mu1", "mu2", "tau1", "tau2", "rho")
  )
  
  if (alpha != 0) {
    H_score <- hyvarinen_score_bnnDPD(alpha, par["mu1"], par["mu2"], Sigma_hat,
                                           yi1, yi2, si1, si2,
                                           jitter = jitter
    )
    H_score <- H_score$score
  } else {
    H_score <- NA
  }
  
  if (compute_variance == TRUE){
    var_ci = add_variance_ci(yi1, yi2, si1, si2, alpha, jitter, par, par_trans_hat)
    cov_theta_nat=var_ci$cov_theta_nat
    est_summary=var_ci$est_summary
    CIs    = data.frame(
      est_summary[est_summary$param %in% c("sensitivity", "specificity"), 
                  c("ci_lower", "ci_upper")],
      row.names=c("sensitivity", "specificity")
    )
    weights_mu <- compute_dpd_weights_mu(
      yi1, yi2, si1, si2,
      mu      = c(m1, m2),
      Sigma   = Sigma_hat,
      alpha   = alpha,
      jitter  = jitter
    )
  } else {
    cov_theta_nat=NA
    est_summary=NA
    CIs=NA
    weights_mu=NA
  }
  ret <- list(
    est_summary = est_summary,
    SeSp = c(sensitivity = pooled_sens, specificity = pooled_spec),
    CIs = CIs,
    par = par,
    par_trans = par_trans_hat,
    Sigma = Sigma_hat,
    cov = cov_theta_nat,
    alpha = alpha,
    H_score = H_score,
    weights_mu = weights_mu,
    opt.obj = opt,
    value = -opt$value, # maximized D (or loglik if alpha~0)
    convergence = opt$convergence,
    message = opt$message,
    call = match.call()
  )
  
  return(ret)
}

dpd_contrib_i <- function(theta, i, Y, Si_list, alpha, jitter = 1e-10) {
  m1 <- theta[1]; m2 <- theta[2]
  a1 <- theta[3]; a2 <- theta[4]; a3 <- theta[5]
  
  mu    <- c(m1, m2)
  Sigma <- make_Sigma(a1, a2, a3)
  
  Vi <- Si_list[[i]] + Sigma
  ld <- logdet_and_solve(Vi, jitter = jitter)
  if (ld$fail) return(-1e12)
  
  diff  <- Y[i, ] - mu
  quad  <- as.numeric(t(diff) %*% ld$inv %*% diff)
  c2pi  <- log(2 * pi)
  logQi <- -c2pi - 0.5 * ld$logdet
  
  if (alpha < 1e-10) {
    # log-likelihood up to constants
    D_i <- -0.5 * (quad + ld$logdet)
  } else {
    log_gi   <- alpha * logQi - 0.5 * alpha * quad
    gi       <- exp(log_gi)
    Qi_alpha <- exp(alpha * logQi)
    D_i      <- gi / alpha - Qi_alpha / ((1 + alpha)^2)
  }
  
  D_i
}

# Variance Computation
compute_sandwich_5d <- function(par_trans_hat, yi1, yi2, si1, si2,
                                alpha, jitter = 1e-10) {
  stopifnot(length(par_trans_hat) == 5)

  n <- length(yi1)
  Y <- cbind(yi1, yi2)
  
  # within-study variance matrices S_i
  Si_list <- lapply(seq_len(n), function(i) {
    diag(c(si1[i]^2, si2[i]^2), nrow = 2)
  })
  
  p <- length(par_trans_hat) # 5
  A_sum <- matrix(0, p, p)
  B_sum <- matrix(0, p, p)
  
  for (i in seq_len(n)) {
    D_i <- function(par) dpd_contrib_i(par, i, Y, Si_list, alpha, jitter)
    
    psi_i <- numDeriv::grad(D_i, par_trans_hat)
    
    Ai    <- numDeriv::hessian(D_i, par_trans_hat)
    
    B_sum <- B_sum + tcrossprod(psi_i)
    A_sum <- A_sum + Ai
  }
  
  A_hat <- A_sum / n
  B_hat <- B_sum / n
  
  A_inv <- tryCatch(solve(A_hat), error = function(e) NULL)
  if (is.null(A_inv)) {
    warning("A_hat is singular; returning NA covariance.")
    cov_theta <- matrix(NA_real_, p, p)
  } else {
    # Sandwich variance: Var(theta_hat) ~= (1/n) A^{-1} B A^{-1}
    cov_theta <- (A_inv %*% B_hat %*% t(A_inv))/n
  }
  
  rownames(cov_theta) <- colnames(cov_theta) <-
    c("mu1", "mu2", "a1", "a2", "a3")
  
  cov_theta
}

transform_cov_to_natural <- function(cov_theta, par_trans_hat) {
  # Delta method
  a1 = par_trans_hat["a1"]; a2 = par_trans_hat["a2"]; a3 = par_trans_hat["a3"]
  rho <- tanh(a3)
  
  J <- diag(5)
  J[3, 3] <- exp(a1)      # d tau1 / d a1
  J[4, 4] <- exp(a2)      # d tau2 / d a2
  J[5, 5] <- 1 - rho^2    # d rho  / d a3 = sech^2(a3)
  
  cov_nat <- J %*% cov_theta %*% t(J)
  rownames(cov_nat) <- colnames(cov_nat) <-
    c("mu1", "mu2", "tau1", "tau2", "rho")
  cov_nat
}

add_variance_ci = function(yi1, yi2, si1, si2, 
                           alpha, jitter, par, par_trans_hat,
                           conf_level = 0.95){
  if (is.null(names(par_trans_hat))) {
    stop("par_trans_hat must be a named vector: c(mu1,mu2,a1,a2,a3)")
  }
  
  # variance/SE estimates
  cov_theta <- compute_sandwich_5d(
    par_trans_hat,
    yi1 = yi1, yi2 = yi2,
    si1 = si1, si2 = si2,
    alpha = alpha,
    jitter = jitter
  )

  cov_theta_nat <- transform_cov_to_natural(cov_theta, par_trans_hat)
  se_nat <- sqrt(diag(cov_theta_nat))
  if (all(is.na(cov_theta_nat))) {
    warning("cov_theta_nat is all NA; SE and CI are set to NA.")
  }

  z <- qnorm(1 - (1 - conf_level)/2)
  lower_ci <- par
  upper_ci <- par

  ## 1) Wald CI on natural scale for mu1, mu2
  idx_mu <- c("mu1", "mu2")
  lower_ci[idx_mu] <- par[idx_mu] - z * se_nat[idx_mu]
  upper_ci[idx_mu] <- par[idx_mu] + z * se_nat[idx_mu]
  
  ## 2) CI for tau1, tau2 via log(tau) = a1,a2 (guarantees tau > 0)
  a1_hat <- par_trans_hat["a1"]
  a2_hat <- par_trans_hat["a2"]
  var_a1 <- cov_theta["a1", "a1"]
  var_a2 <- cov_theta["a2", "a2"]
  
  if (is.na(var_a1) || var_a1 < 0) {
    warning("Variance for a1 is invalid; tau1 CI set to NA.")
    lower_ci["tau1"] <- NA_real_
    upper_ci["tau1"] <- NA_real_
  } else {
    se_a1 <- sqrt(var_a1)
    lower_a1 <- a1_hat - z * se_a1
    upper_a1 <- a1_hat + z * se_a1
    lower_ci["tau1"] <- exp(lower_a1)
    upper_ci["tau1"] <- exp(upper_a1)
  }
  
  if (is.na(var_a2) || var_a2 < 0) {
    warning("Variance for a2 is invalid; tau2 CI set to NA.")
    lower_ci["tau2"] <- NA_real_
    upper_ci["tau2"] <- NA_real_
  } else {
    se_a2 <- sqrt(var_a2)
    lower_a2 <- a2_hat - z * se_a2
    upper_a2 <- a2_hat + z * se_a2
    lower_ci["tau2"] <- exp(lower_a2)
    upper_ci["tau2"] <- exp(upper_a2)
  }

  ## 2) CI for rho via a3-scale (guarantees [-1, 1])
  a3_hat <- par_trans_hat["a3"]
  var_a3 <- cov_theta["a3", "a3"]
  if (is.na(var_a3) || var_a3 < 0) {
    warning("Variance for a3 is invalid; rho CI set to NA.")
    lower_ci["rho"] <- NA_real_
    upper_ci["rho"] <- NA_real_
    se_rho <- NA_real_
  } else {
    se_a3 <- sqrt(var_a3)
    # Make Wald-type CI on a3 scale and transform endpoints back to rho = tanh(a3)
    lower_a3 <- a3_hat - z * se_a3
    upper_a3 <- a3_hat + z * se_a3
    lower_ci["rho"] <- tanh(lower_a3)
    upper_ci["rho"] <- tanh(upper_a3)
  }

  se_par <- c(
    se_nat["mu1"],
    se_nat["mu2"],
    se_nat["tau1"],
    se_nat["tau2"],
    se_nat["rho"]
  )
  names(se_par) <- c("mu1", "mu2", "tau1", "tau2", "rho")

  # confidence intervals for sensitivity/specificity
  lower_ci_sens <- plogis(lower_ci["mu1"])
  lower_ci_spec <- plogis(lower_ci["mu2"])
  upper_ci_sens <- plogis(upper_ci["mu1"])
  upper_ci_spec <- plogis(upper_ci["mu2"])
  # SE for sensitivity/specificity (using delta method)
  p_sens=plogis(par["mu1"])
  p_spec=plogis(par["mu2"])
  se_sens <- p_sens * (1 - p_sens) * se_par["mu1"]
  se_spec <- p_spec * (1 - p_spec) * se_par["mu2"]
  
  # ============================================================
  # HKSJ adjustment for (mu1, mu2) only
  #  - hksj: t(df=2N-2) only
  # ============================================================
  n <- length(yi1)
  DF <- 2*n - 2
  tcrit <- qt(1 - (1 - conf_level) / 2, df = DF)
  
  # Build Sigma_hat from natural parameters (tau1, tau2, rho)
  # (assumes tau1,tau2 are SDs, not variances)
  tau1_hat <- par["tau1"]
  tau2_hat <- par["tau2"]
  rho_hat  <- par["rho"]
  
  Sigma_hat <- matrix(c(tau1_hat^2, rho_hat * tau1_hat * tau2_hat,
                        rho_hat * tau1_hat * tau2_hat, tau2_hat^2),
                      2, 2, byrow = TRUE)
  
  mu_hat <- unname(par[idx_mu])
  
  # HKSJ CIs on mu-scale
  lower_ci_hksj_mu <- mu_hat - tcrit * se_nat[idx_mu]
  upper_ci_hksj_mu <- mu_hat + tcrit * se_nat[idx_mu]
  
  names(lower_ci_hksj_mu) <- names(upper_ci_hksj_mu) <- idx_mu
  
  # Also provide HKSJ CIs on probability scale (Se/Sp) by transforming mu endpoints
  lower_ci_hksj_sens <- plogis(lower_ci_hksj_mu["mu1"])
  upper_ci_hksj_sens <- plogis(upper_ci_hksj_mu["mu1"])
  lower_ci_hksj_spec <- plogis(lower_ci_hksj_mu["mu2"])
  upper_ci_hksj_spec <- plogis(upper_ci_hksj_mu["mu2"])
  
  # Build summary table
  est_summary <- data.frame(
    param    = c(names(par), "sensitivity", "specificity"),
    estimate = c(par, p_sens, p_spec),
    se       = c(se_par, se_sens, se_spec),
    ci_lower = c(lower_ci, lower_ci_sens, lower_ci_spec),
    ci_upper = c(upper_ci, upper_ci_sens, upper_ci_spec),
    ci_lower_hksj = NA_real_,
    ci_upper_hksj = NA_real_,
    row.names = NULL
  )
  
  # Fill HKSJ columns for mu1/mu2 rows
  row_mu1 <- which(est_summary$param == "mu1")
  row_mu2 <- which(est_summary$param == "mu2")
  est_summary$ci_lower_hksj[row_mu1] <- lower_ci_hksj_mu["mu1"]
  est_summary$ci_upper_hksj[row_mu1] <- upper_ci_hksj_mu["mu1"]
  
  est_summary$ci_lower_hksj[row_mu2] <- lower_ci_hksj_mu["mu2"]
  est_summary$ci_upper_hksj[row_mu2] <- upper_ci_hksj_mu["mu2"]
  
  # Fill HKSJ columns for sensitivity/specificity rows (prob-scale)
  row_sens <- which(est_summary$param == "sensitivity")
  row_spec <- which(est_summary$param == "specificity")
  
  est_summary$ci_lower_hksj[row_sens] <- lower_ci_hksj_sens
  est_summary$ci_upper_hksj[row_sens] <- upper_ci_hksj_sens
  
  est_summary$ci_lower_hksj[row_spec] <- lower_ci_hksj_spec
  est_summary$ci_upper_hksj[row_spec] <- upper_ci_hksj_spec
  
  return(list(est_summary = est_summary,
              cov_theta_nat = cov_theta_nat) 
  )
}

# study-specific contribution measures computation
compute_dpd_weights_mu <- function(yi1, yi2, si1, si2, mu, Sigma, alpha, jitter = 1e-10) {
  Y <- cbind(yi1, yi2)
  n <- nrow(Y)
  
  # Within-study variance matrices are diagonal
  Si_list <- lapply(seq_len(n), function(i) diag(c(si1[i]^2, si2[i]^2), nrow = 2))
  
  p <- length(mu)
  c2pi <- log(2 * pi)
  
  # store C_i = g_i W_i and accumulate C = sum C_i
  C_list <- vector("list", n)
  C_sum <- matrix(0, nrow = p, ncol = p)
  
  # diagnostics: g_i
  g <- rep(NA_real_, n)
  
  for (i in seq_len(n)) {
    Vi <- Si_list[[i]] + Sigma
    ld <- logdet_and_solve(Vi, jitter = jitter)
    if (ld$fail) {
      stop(sprintf("Non-PD Vi at study %d in weight calculation.", i))
    }
    
    r <- Y[i, ] - mu
    Wi <- ld$inv
    
    quad <- as.numeric(t(r) %*% Wi %*% r)
    logQi <- -c2pi - 0.5 * ld$logdet
    
    # g_i
    gi <- exp(alpha * logQi - 0.5 * alpha * quad)
    g[i] <- gi
    
    Ci <- gi * Wi
    C_list[[i]] <- Ci
    C_sum <- C_sum + Ci
  }
  
  C_sum_inv <- logdet_and_solve(C_sum, jitter = jitter)
  if (C_sum_inv$fail) {
    warning("Matrix B is singular in effective-weight computation; returning NA outputs.")
    U_list <- replicate(n, matrix(NA_real_, nrow = p, ncol = p), simplify = FALSE)
    U_mu1 <- U_mu2 <- rep(NA_real_, n)
  } else {
    # weight matrices
    U_list <- lapply(C_list, function(Ci) C_sum_inv$inv %*% Ci)
    
    u_avg <- sapply(U_list, function(Ui) 
      (abs(Ui[1, 1])+abs(Ui[2, 2])+abs(Ui[1, 2])+abs(Ui[2,1]))
      )
    u_avg <- u_avg/sum(u_avg)
    u_mu1 <- sapply(U_list, function(Ui) abs(Ui[1, 1]))
    u_mu2 <- sapply(U_list, function(Ui) abs(Ui[2, 2]))
    u_mu1 <- u_mu1/sum(u_mu1)
    u_mu2 <- u_mu2/sum(u_mu2)
  }
  
  list(
    U_list = U_list,
    u_avg  = u_avg,
    u_mu1  = u_mu1,
    u_mu2  = u_mu2,
    g = g
  )
}

# Hyvarinen score for DPD based on a fitted robust BNN model
# Computes per-study and aggregate scores
hyvarinen_score_bnnDPD <- function(alpha, mu1, mu2, Sigma,
                                        yi1, yi2, si1, si2, jitter = 1e-10) {
  stopifnot(
    length(yi1) == length(yi2),
    length(yi1) == length(si1),
    length(yi1) == length(si2)
  )

  n <- length(yi1)
  Y <- cbind(yi1, yi2)

  mu <- c(mu1, mu2)

  # Prepare within-study variances S_i
  Si_list <- lapply(seq_len(n), function(i) diag(c(si1[i]^2, si2[i]^2), nrow = 2))

  # Constants
  p <- 2L
  c2pi <- log(2 * pi)

  H <- numeric(n)

  for (i in 1:n) {
    Vi <- Si_list[[i]] + Sigma
    ld <- logdet_and_solve(Vi, jitter = jitter)
    if (ld$fail) stop(sprintf("Non-PD (after jitter) Vi at study %d", i))

    ri <- as.numeric(Y[i, ] - mu)

    # u = W r via triangular solves using the Cholesky factor
    # chol is upper-triangular R s.t. t(R) %*% R = V
    z <- forwardsolve(t(ld$chol), ri)
    u <- backsolve(ld$chol, z)
    # Quadratic form Q = r' W r
    Q <- sum(ri * u)

    # trace(W_i)
    trW <- sum(diag(ld$inv))
    # r' W_i W_i r = || W r ||^2
    rW2r <- sum(u * u)

    # Density f_i(y_i) under N(mu, V_i)
    logfi <- -0.5 * (p * c2pi + ld$logdet + Q)
    log_f_a  <- alpha * logfi
    log_f_2a <- 2 * alpha * logfi
    f_a  <- exp(log_f_a)
    f_2a <- exp(log_f_2a)

    H[i] <- f_a * (alpha * rW2r - trW) + 0.5 * f_2a * rW2r
  }

  list(
    score    = sum(H),
    scores_i = H
  )
}

# Optimize alpha by minimizing the Hyvarinen score (wrapper around fit_bnnDPD_core)
#
# Strategy (default): grid search over alpha in (alpha_min, alpha_max),
# using warm-starts from the previous alpha's unconstrained parameters for speed and stability.
# Returns the best fit and the score path.
bnnDPD_grid_search <- function(yi1, yi2, si1, si2,
                               gamma = 0.01,
                               alpha_lim = c(0.01, 0.33),
                               init.par = NULL,
                               seed = NULL,
                               control = list(maxit = 1000, reltol = 1e-4),
                               jitter = 1e-10, verbose = TRUE) {
  stopifnot(gamma > 0, gamma < 1)
  stopifnot(length(alpha_lim) == 2, 
            alpha_lim[1] < alpha_lim[2], 
            alpha_lim[1] >= 0)
  
  alphas <- seq(alpha_lim[1], alpha_lim[2], by = gamma)

  fits <- vector("list", length(alphas))
  scores <- rep(NA_real_, length(alphas))

  # We'll warm-start: pass the previous fit's unconstrained params (opt.obj$par)
  init_par_curr <- init.par

  if (verbose) {
    cat(sprintf("Searching alpha on grid [%g, %g] with %d points\n",
                min(alphas), max(alphas), length(alphas)))
  }

  for (num in seq_along(alphas)) {
    alpha <- alphas[num]
    
      fit <- fit_bnnDPD_core(yi1, yi2, si1, si2,
                                 alpha = alpha,
                                 init.par = init_par_curr,
                                 control = control,
                                 jitter = jitter,
                                 compute_variance=FALSE,
                                 verbose = FALSE)

    fits[[num]] <- fit
    scores[num] <- fit$H_score

    # warm start for next alpha
    if (!is.null(fit$opt.obj) && !is.null(fit$opt.obj$par)) {
      init_par_curr <- c(fit$par["mu1"],fit$par["mu2"],fit$opt.obj$par)
    }
  }

  # pick alpha minimizing Hyvarinen score
  finite_scores <- ifelse(is.finite(scores), scores, Inf)
  best_idx <- which.min(finite_scores)
  best_fit <- fits[[best_idx]]
  
  if (verbose) {
    cat(sprintf("Best alpha = %.4f with H-score = %.6f\n",
                alphas[best_idx], scores[best_idx]))
  }
  
  # add estimates for variance, ci,and study weights 
  var_ci = add_variance_ci(yi1, yi2, si1, si2, 
                           alpha=best_fit$alpha, jitter=jitter,
                           par = best_fit$par,
                           par_trans_hat = best_fit$par_trans)
  
  weights_mu <- compute_dpd_weights_mu(
    yi1, yi2, si1, si2,
    mu      = c(best_fit$par["mu1"], best_fit$par["mu2"]),
    Sigma   = best_fit$Sigma,
    alpha   = best_fit$alpha,
    jitter  = jitter
  )
  
  est_summary=var_ci$est_summary
  best_fit$cov_theta_nat=var_ci$cov_theta_nat
  best_fit$est_summary=est_summary
  best_fit$weights_mu=weights_mu
  
  tbl <- do.call(rbind, lapply(fits, function(fit) {
    c(alpha = fit$alpha, unclass(fit$par), unclass(fit$SeSp), H_score = fit$H_score)
  }))
  tbl <- as.data.frame(tbl, check.names = FALSE)
  rownames(tbl) <- NULL
  
  list(
    SeSp   = best_fit$SeSp,
    Sigma  = best_fit$Sigma,
    CIs    = data.frame(
      est_summary[est_summary$param %in% c("sensitivity", "specificity"), 
                  c("ci_lower", "ci_upper")],
      row.names=c("sensitivity", "specificity")
      ),
    est_summary = best_fit$est_summary,
    par = best_fit$par,
    weights_mu = best_fit$weights_mu,
    conv = best_fit$convergence,
    best_fit = best_fit,
    best_alpha = alphas[best_idx],
    scores = scores,
    alphas = alphas,
    fits = fits,
    est_table = tbl
  )
}

# fit BNN-DPD with optimize alpha by minimizing the Hyvarinen score
# initial values are set to estimates from the standard BNN model
fit_bnnDPD_grid_search <- function(dat, 
                               gamma=0.01, alpha_lim=c(0.01, 0.33),
                               verbose=FALSE){
  df <- make_logit_data(TP=dat$TP, FN=dat$FN, FP=dat$FP, TN=dat$TN, 
                        correction.control = "single")
  
  # --- initial values from standard BNN (alpha = 0) ---
  res_bnn = fit_bnn(Data = dat, method = "ml")
  
  check_mu  <- any(!is.finite(res_bnn$SeSp.logit))
  check_Psi <- any(!is.finite(diag(res_bnn$Psi)))
  if (check_mu) {
    warning("fit_bnn returned non-finite SeSp.logit; default initial values are used.")
  }
  if (check_Psi) {
    warning("fit_bnn returned non-finite diagonal elements of Psi; default initial values are used.")
  }
  
  if (check_mu || check_Psi){
    init.par_trans <- NULL
  } else {
    MLE.par_nat <- c(
      res_bnn$SeSp.logit,                                    # mu1, mu2
      sqrt(res_bnn$Psi[1, 1]),                               # tau1
      sqrt(res_bnn$Psi[2, 2]),                               # tau2
      res_bnn$Psi[1, 2] / sqrt(res_bnn$Psi[1, 1] * res_bnn$Psi[2, 2]) # rho
    )
    names(MLE.par_nat) <- c("mu1", "mu2", "tau1", "tau2", "rho")
    
    tau1_safe <- max(MLE.par_nat["tau1"], 1e-4)
    tau2_safe <- max(MLE.par_nat["tau2"], 1e-4)
    rho_clipped <- min(max(MLE.par_nat["rho"], -0.999), 0.999)
    
    init.par_trans <- c(
      mu1 = MLE.par_nat["mu1"],
      mu2 = MLE.par_nat["mu2"],
      a1  = log(tau1_safe),
      a2  = log(tau2_safe),
      a3  = atanh(rho_clipped)
    )
  }
  
  fit_DPD = bnnDPD_grid_search(df$yi1, df$yi2, df$si1, df$si2, 
                               init.par = init.par_trans,
                               gamma=gamma, alpha_lim=alpha_lim,
                               verbose = verbose)
  return (fit_DPD)
}

plot_dpd_contribution_rate <- function(res_obj,
                                       ylim_g = NULL,
                                       ylim_u = c(0, 1),
                                       facet_labels_u = c("Sensitivity", "Specificity"),
                                       percent_accuracy = 1) {
  
  stopifnot(!is.null(res_obj$weights_mu))
  stopifnot(all(c("g", "u_mu1", "u_mu2") %in% names(res_obj$weights_mu)))
  
  g     <- res_obj$weights_mu$g
  u_mu1 <- res_obj$weights_mu$u_mu1
  u_mu2 <- res_obj$weights_mu$u_mu2
  
  k <- length(g)
  if (length(u_mu1) != k || length(u_mu2) != k) stop("Lengths of g, u_mu1, u_mu2 must match.")
  
  theme_paper <- ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 11),
      axis.text  = ggplot2::element_text(size = 11),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
  
  # ---- Left: g plot (raw scale) ----
  df_g <- tibble::tibble(
    Study = seq_len(k),
    value = g
  )
  
  p_g <- ggplot2::ggplot(df_g, ggplot2::aes(x = Study, y = value)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.55) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(k),
      limits = c(1, k),
      expand = ggplot2::expansion(mult = c(0.02, 0.03))
    ) +
    ggplot2::coord_cartesian(ylim = ylim_g) +
    ggplot2::labs(x = "Study", 
                  y = "DPD weight")+
    theme_paper
  
  # ---- Right: u_mu1/u_mu2 facet plot (percent scale) ----
  df_u <- tibble::tibble(
    Study = seq_len(k),
    u_mu1 = u_mu1,
    u_mu2 = u_mu2
  ) |>
    tidyr::pivot_longer(cols = c(u_mu1, u_mu2),
                        names_to = "component",
                        values_to = "value") |>
    dplyr::mutate(
      component = factor(
        component,
        levels = c("u_mu1", "u_mu2"),
        labels = facet_labels_u
      )
    )
  
  p_u <- ggplot2::ggplot(df_u, ggplot2::aes(x = Study, y = value)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.55) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::facet_wrap(~ component, nrow = 1, strip.position = "top") +
    ggplot2::scale_x_continuous(
      breaks = seq_len(k),
      limits = c(1, k),
      expand = ggplot2::expansion(mult = c(0.02, 0.03))
    ) +
    ggplot2::scale_y_continuous(
      limits = ylim_u,
      labels = scales::percent_format(accuracy = percent_accuracy)
    ) +
    ggplot2::labs(x = "Study", y = "Contribution rate (%)") +
    ggplot2::theme(
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    ) +
    theme_paper
  
  p_panel <- patchwork::wrap_plots(
    patchwork::free(p_g, side = "b"),
    p_u,
    nrow = 1,
    widths = c(1, 2)
  )
  
  p_all =list(p_g=p_g, p_u=p_u, p_panel=p_panel)
  return(p_all)
}

# The BNN model of Reitsma et al. (2005)
fit_bnn <- function(Data, method=c("reml", "ml")){
  # "Data" is a list of data frames consisting of the 2X2 DTA table from several studies.
  # The column names of "Data" must include the frequencies "TP", "FN", "FP" and "TN".
  
  k <- nrow(Data)
  
  # Compute the observed logit-transformed Se and Sp and the within-study covariance matrices
  Y <- comp.y(Data = Data)
  Psi <- comp.Psi(Data = Data)
  
  correction <- 0.5
  for(i in 1:k){
    if(Data$TP[i]==0|Data$FP[i]==0|Data$FN[i]==0|Data$TN[i]==0){
      Data$TP[i] = Data$TP[i] + correction
      Data$FP[i] = Data$FP[i] + correction
      Data$FN[i] = Data$FN[i] + correction
      Data$TN[i] = Data$TN[i] + correction
    }
  }
  
  TP <- Data$TP; FP <- Data$FP; FN <- Data$FN; TN <- Data$TN  
  n1 <- TP + FN; n2 <- FP + TN; npat <- n1 + n2
  
  cov.sesp.logit <- rep(0, k)
  trsesp.logit <- matrix(0, nrow = k, ncol = 2)
  trsesp.logit <- cbind(tsens = escalc(measure="PLO", xi=TP, ni=n1)[,1],
                        tfpr = escalc(measure="PLO", xi=TN, ni=n2)[,1])
  
  var.se.logit = var.sp.logit <- rep(0, k)
  var.se.logit <- escalc(measure="PLO", xi=TP, ni=n1)[,2]
  var.sp.logit <- escalc(measure="PLO", xi=TN, ni=n2)[,2]
  
  Data.mvmeta.logit <- matrix(0, nrow = k, ncol=6)
  Data.mvmeta.logit <- data.frame(pat.num=npat, trse=trsesp.logit[,1], trsp=trsesp.logit[,2],
                                  varse=var.se.logit, covsefpr=cov.sesp.logit, varsp=var.sp.logit)
  
  fit.logit <- list()
  mvmeta.logit <- mvmeta(cbind(trse, trsp)~1, S=Data.mvmeta.logit[, 4:6], method = method, data = Data.mvmeta.logit, control=list(maxiter=1000))
  
  summary.logit <- summary(mvmeta.logit)$coefficients
  SeSp <- c(ilogit(summary.logit[1,1]), ilogit(summary.logit[2,1]))
  SeCI <- as.numeric(ilogit(summary.logit[1,5:6]))
  SpCI <- as.numeric(ilogit(summary.logit[2,5:6]))
  
  loglik.bnn <- mvmeta.logit$logLik
  
  fit.logit$mvmeta.logit <- summary(mvmeta.logit)
  fit.logit$SeSp <- SeSp
  fit.logit$SeSp.logit <- logit(SeSp)
  fit.logit$Psi <- mvmeta.logit$Psi
  fit.logit$vcov <- mvmeta.logit$vcov
  fit.logit$CIs <- c("Se.lb"=SeCI[1],"Se.ub"=SeCI[2], "Sp.lb"=SpCI[1],"Sp.ub"=SpCI[2])
  fit.logit$Stats <- c(logLik=loglik.bnn, AIC=-2*loglik.bnn + 2*5, BIC=-2*loglik.bnn + 5*log(k*2))
  
  return(fit.logit)
  
}

# 2x2 -> logit & SE transform
make_logit_data <- function(TP, FN, FP, TN,
                            correction = 0.5,
                            correction.control = c("all", "single", "none")) {
  # match correction option
  correction.control <- match.arg(correction.control)
  
  stopifnot(
    length(TP) == length(FN),
    length(TP) == length(FP),
    length(TP) == length(TN)
  )
  
  TP <- as.numeric(TP)
  FN <- as.numeric(FN)
  FP <- as.numeric(FP)
  TN <- as.numeric(TN)
  
  # continuity correction
  if (correction.control == "all") {
    # if any cell is zero in any study, add correction to all cells of all studies
    has_zero <- (TP == 0) | (FN == 0) | (FP == 0) | (TN == 0)
    if (any(has_zero, na.rm = TRUE)) {
      TP <- TP + correction
      FN <- FN + correction
      FP <- FP + correction
      TN <- TN + correction
    }
    
  } else if (correction.control == "single") {
    # add correction only to studies that have at least one zero cell
    zero_study <- (TP == 0) | (FN == 0) | (FP == 0) | (TN == 0)
    zero_study[is.na(zero_study)] <- FALSE  # do not treat NA as zero
    
    if (any(zero_study)) {
      TP[zero_study] <- TP[zero_study] + correction
      FN[zero_study] <- FN[zero_study] + correction
      FP[zero_study] <- FP[zero_study] + correction
      TN[zero_study] <- TN[zero_study] + correction
    }
    # "none" case: no correction
  }
  
  # probabilities
  sens <- TP / (TP + FN)
  spec <- TN / (FP + TN)
  
  # logit transform
  yi1 <- qlogis(sens)  # logit(sensitivity)
  yi2 <- qlogis(spec)  # logit(specificity)
  
  # delta method
  vi1 <- 1 / TP + 1 / FN
  vi2 <- 1 / TN + 1 / FP
  
  data.frame(
    yi1 = yi1,
    yi2 = yi2,
    si1 = sqrt(vi1),
    si2 = sqrt(vi2)
  )
}

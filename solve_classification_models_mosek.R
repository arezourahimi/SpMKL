library(Rmosek)

solve_classification_svm <- function(K, y, C, epsilon) {
  N <- length(y)
  yyK <- (y %*% t(y)) * K
  
  problem <- list()
  problem$sense <- "min"
  problem$c <- rep(-1, N)
  problem$A <- Matrix(y, nrow = 1, byrow = TRUE, sparse = TRUE)
  problem$bc <- rbind(blc = 0, buc = 0) 
  
  bux <- rep(C, N)
  
  problem$bx <- rbind(blx = rep(0, N), bux = bux)
  
  I <- matrix(1:N, N, N, byrow = FALSE)
  J <- matrix(1:N, N, N, byrow = TRUE)
  problem$qobj <- list(i = I[lower.tri(I, diag = TRUE)],
                       j = J[lower.tri(J, diag = TRUE)],
                       v = yyK[lower.tri(yyK, diag = TRUE)])
  
  problem$iparam <- list(INTPNT_MULTI_THREAD=0, NUM_THREADS=1)
  
  
  opts <- list()
  opts$verbose <- 0
  start_time <- Sys.time()
  result <- mosek(problem, opts)
  end_time <- Sys.time()
  elapsed_time <- as.numeric(end_time-start_time, units="secs")
  
  alpha <- result$sol$itr$xx[1:N]
  alpha_original <- alpha
  alpha[alpha < +C * epsilon] <- 0
  alpha[alpha > +C * (1 - epsilon)] <- +C
  
  objective <- sum(alpha) - 0.5 * (t(alpha) %*% yyK) %*% (alpha)
  objective <- objective * (objective >= 0)
  
  objective_original <- sum(alpha_original) - 0.5 * (t(alpha_original) %*% yyK) %*% (alpha_original)
  objective_original <- objective_original * (objective_original >= 0)
  
  support_indices <- which(alpha != 0)
  active_indices <- which(alpha != 0 & alpha < C)
  
  if (length(active_indices) > 0) {
    b <- mean(y[active_indices] * (1 - yyK[active_indices, support_indices] %*% alpha[support_indices]))
  } else {
    b <- 0
  }
  
  model <- list(alpha = alpha * y, b = b, objective = objective, alpha_original = alpha_original*y, objective_original = objective_original, Time = elapsed_time)
  
  return(model)
}

solve_master_problem <- function(obj_coef, constraints_matrix, rhs, senses, lb, ub, vtype, P, contains_g = TRUE, LP_relaxation = FALSE) {
  
  solve_as_lp <- LP_relaxation | !contains_g
  
  if(solve_as_lp)
    vtype <- rep("C", length(vtype))
  
  blc <- rhs
  buc <- rhs
  blc[which(senses=="L")] <- -Inf
  buc[which(senses=="G")] <- Inf
  
  problem <- list()
  problem$sense <- "min"
  problem$c <- obj_coef
  problem$A <- Matrix(constraints_matrix)
  problem$bc <- rbind(blc = blc, buc = buc) 
  problem$bx <- rbind(blx = lb, bux = ub)
  
  if(!solve_as_lp)
    problem$intsub <- which(vtype=="B")
  
  problem$iparam <- list(INTPNT_MULTI_THREAD=0, NUM_THREADS=1)
  
  opts <- list()
  opts$verbose <- 0
  
  start_time <- Sys.time()
  result <- mosek(problem, opts)
  end_time <- Sys.time()
  elapsed_time <- as.numeric(end_time-start_time, units="secs")
  
  if(solve_as_lp)
  {  
    sln <- result$sol$itr$xx
  }else
  {
    sln <- result$sol$int$xx
  }

  Gamma <- sln[1]
  eta <- sln[2:(P+1)]
  
  objective <- Gamma
  
  output <- list(Gamma = Gamma, eta = eta, objective = objective, Time = elapsed_time)
  
  if(contains_g)
    output$g <- sln[(P+2):(2*P+1)]
  
  return(output)
}
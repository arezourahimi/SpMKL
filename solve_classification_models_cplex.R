library(Rcplex)

solve_classification_svm <- function(K, y, C, epsilon) {
  N <- length(y)
  yyK <- (y %*% t(y)) * K
  
  opts <- list()
  opts$trace <- 0
  opts$maxcalls <- 41 * 200
  
  start_time <- Sys.time()
  result <- Rcplex(cvec = rep(-1, N), 
                   Amat = matrix(y, nrow = 1, byrow = TRUE), 
                   bvec = 0, 
                   Qmat = yyK, 
                   lb = rep(0, N),
                   ub = rep(C, N),
                   control = opts,
                   objsense = "min",
                   sense = "E")
  end_time <- Sys.time()
  elapsed_time <- as.numeric(end_time-start_time, units="secs")
  
  alpha <- result$xopt[1:N]
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
  
  opts <- list()
  opts$trace <- 0
  opts$maxcalls <- 41 * 200
  
  start_time <- Sys.time()
  result <- Rcplex(cvec = obj_coef, 
                   Amat = constraints_matrix, 
                   bvec = rhs, 
                   lb = lb,
                   ub = ub,
                   vtype = vtype,
                   control = opts,
                   objsense = "min",
                   sense = senses)
  end_time <- Sys.time()
  elapsed_time <- as.numeric(end_time-start_time, units="secs")
  
  sln <- result$xopt
  Gamma <- sln[1]
  eta <- sln[2:(P+1)]
  
  objective <- Gamma
  
  output <- list(Gamma = Gamma, eta = eta, objective = objective, Time = elapsed_time)
  
  if(contains_g)
    output$g <- sln[(P+2):(2*P+1)]
  
  return(output)
}
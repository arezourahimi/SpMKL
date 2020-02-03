spmkl_step1_phase2_cutting_plane <- function(Km, y, parameters, initial_info = NULL) {
  # parameters should contain: G
  
  P <- dim(Km)[3]
  
  terminate <- FALSE
  LB <- -Inf
  
  #incumbent solution
  incumbent_solution <- list(UB = Inf, eta = NULL, g = NULL, alpha = NULL, b = NULL)
  
  iteration <- 1
  
  #define constraint matrix here, and add the equality constraint. Also define the objective function
  #order of variables: Gamma, eta_1,...,eta_P, g_1,...,g_P
  NVars <- 1 + 2*P
  master_obj_coef <- numeric(NVars)
  master_obj_coef[1] <- 1
  
  master_constraint_matrix <- matrix(0, nrow = 2 + P, ncol = NVars)
  master_rhs <- rep(0,2+P)
  master_senses <- rep("L",2+P)
  master_lb <- rep(0, NVars)
  master_ub <- rep(1, NVars)
  master_vtype <- rep("B", NVars)
  
  master_vtype[1:(P+1)] <- "C"
  master_ub[1] <- Inf
  
  # sum eta = 1
  master_constraint_matrix[1, 2:(P+1)] <- 1
  master_rhs[1] <- 1
  master_senses[1] <- "E"
  
  # sum g <= G
  master_constraint_matrix[2, (2+P):(1+2*P)] <- 1
  master_rhs[2] <- parameters$G
  
  # eta_m <= g_m
  for(m in 1:P)
    master_constraint_matrix[2+m, c(m+1, m+1+P)] <- c(1,-1)
  
  #we store the newly generated cuts, in case they are needed in the future  
  new_svm_constraints <- matrix(0, nrow = 0, ncol = NVars)
  new_svm_rhs <- c()
  
  J <- 0
  Gamma <- 0 #approximator of J
  
  alpha <- NULL
  b <- 0
  eta <- NULL
  g <- NULL
  
  initial_cuts_added <- 0
  
  if(is.null(initial_info) == FALSE)
  {
    if(!is.null(initial_info$svm_cuts))
    {
      initial_cuts_added <- length(initial_info$svm_cuts$rhs)
      initial_cuts <- matrix(0, nrow = initial_cuts_added, ncol = NVars)
      initial_cuts[,1:(P+1)] <- initial_info$svm_cuts$matrix[,1:(P+1)]
      master_constraint_matrix <- rbind(master_constraint_matrix, initial_cuts)
      master_rhs <- c(master_rhs, initial_info$svm_cuts$rhs)
      master_senses <- c(master_senses, rep("G", initial_cuts_added))
      
      print(sprintf("Added %d initial SVM cuts.", initial_cuts_added))
    }
    
    if(!is.null(initial_info$solution))
    {
      eta <- initial_info$solution$eta
      if(parameters$G > 0)
        eta <- sparsify_eta(eta, parameters$G)
      
      g <- 1*(eta>0)
      b <- initial_info$solution$b
      alpha <- initial_info$solution$alpha
      J <- initial_info$solution$J
      
      incumbent_solution <- list(UB = J, eta = eta, g = g, alpha = alpha, b = b)
    }
  }
  
  objectives <- c()
  
  info_cols <- c("obj", "UB", "LB", "Gap", "Gamma", "J", "Time(M)", "Time(SVM)")
  iteration_info <- matrix(0, 0, length(info_cols))
  colnames(iteration_info) <- info_cols
  
  total_time <- 0
  
  while(terminate == FALSE)
  {
    if(iteration > 1 | is.null(eta) | initial_cuts_added > 0)
    {
      #SOLVE MP and obtain Gamma, eta and g
      master_result <- solve_master_problem(master_obj_coef, master_constraint_matrix, master_rhs, master_senses, 
                                                   master_lb, master_ub, master_vtype, P)
      Gamma <- master_result$Gamma
      eta <- master_result$eta
      g <- master_result$g
      master_time <- master_result$Time
    }else
    {
      Gamma <- 0
      master_time <- 0
    }
    
    Keta <- calculate_Keta(Km, eta)
    svm_result <- solve_classification_svm(Keta, y, parameters$C, parameters$epsilon)
    
    alpha <- svm_result$alpha #regularized alpha (may be numerically infeasible)
    b <- svm_result$b
    
    alpha_cut <- svm_result$alpha_original #numerically feasible
    J <- svm_result$objective_original
    
    #model$alpha_original is already multiplied by y
    alpha_sum <- t(alpha_cut) %*% y
    kappa <- calucalte_kappa(alpha_cut, Km, P)
    
    #add constraint Gamma_t + Eta_t*kappa >= alpha_sum
    constraint_alpha <- numeric(NVars)
    constraint_alpha[1] <- 1
    constraint_alpha[2:(P+1)] <- kappa
    
    master_constraint_matrix <- rbind(master_constraint_matrix, constraint_alpha)
    master_rhs <- c(master_rhs, alpha_sum)
    master_senses <- c(master_senses, "G")  
    
    new_svm_constraints <- rbind(new_svm_constraints, constraint_alpha)
    new_svm_rhs <- c(new_svm_rhs, alpha_sum)
    
    svm_time <- svm_result$Time
    
    total_time <- total_time + master_time + svm_time
    
    obj <- J
    if(obj < incumbent_solution$UB)
    {
      UB <- obj
      incumbent_solution = list(UB=obj, eta = eta, g = g, alpha = alpha, b = b)
    }
    
      LB <- Gamma
    
    Gap <- (incumbent_solution$UB-LB)/max(abs(LB),parameters$epsilon)
    
    objectives <- c(objectives, obj)
    
    info <- c(obj, incumbent_solution$UB, LB, Gap, Gamma, J, master_time, svm_time)
    iteration_info <- rbind(iteration_info, info)
    
    if(iteration %% 5 == 1)
      print(sprintf("%d %f %f %f %f%% %f %f %f %f",iteration, obj, incumbent_solution$UB, LB, Gap*100, Gamma, J, master_time, svm_time))
    
    if(Gap <= parameters$optimality_gap | iteration >= parameters$iteration_count)
    {
      terminate <- TRUE
    }
    
    if(terminate & iteration %% 5 != 1)
      print(sprintf("%d %f %f %f %f%% %f %f %f %f",iteration, obj, incumbent_solution$UB, LB, Gap*100, Gamma, J, master_time, svm_time))
    
    iteration <- iteration + 1
  }
  
  reta <- polish_eta(incumbent_solution$eta, parameters$epsilon)
  
  state <- list(alpha = incumbent_solution$alpha, b = incumbent_solution$b, eta = reta, g = incumbent_solution$g, UB = incumbent_solution$UB,
                objectives = objectives, iteration_info = iteration_info, parameters = parameters, total_time = total_time)
  
  new_cuts_info = list(SVM_cuts = list(matrix = new_svm_constraints, rhs = new_svm_rhs))
  
  output <- list(state = state, new_cuts_info = new_cuts_info)
  
  return(output)
}

spmkl_step1_phase13_cutting_plane <- function(Km, y, parameters, initial_info = NULL) {
  # parameters should contain: G
  
  G <- 0
  if(parameters$sparsify == TRUE)
    G <- parameters$G
  
  P <- dim(Km)[3]
  
  terminate <- FALSE
  LB <- -Inf
  
  #incumbent solution
  incumbent_solution <- list(UB = Inf, eta = NULL, alpha = NULL, b = NULL)
  
  iteration <- 1
  
  #define constraint matrix here, and add the equality constraint. Also define the objective function
  #order of variables: Gamma, eta_1,...,eta_P
  NVars <- 1 + P
  master_obj_coef <- numeric(NVars)
  master_obj_coef[1] <- 1
  
  master_constraint_matrix <- matrix(0, nrow = 1, ncol = NVars)
  master_rhs <- rep(1,1)
  master_senses <- rep("E",1)
  master_lb <- rep(0, NVars)
  master_ub <- rep(1, NVars)
  master_vtype <- rep("C", NVars)
  master_ub[1] <- Inf
  
  # sum eta = 1
  master_constraint_matrix[1, 2:(P+1)] <- 1
  master_rhs[1] <- 1
  master_senses[1] <- "E"
  
  #we store the newly generated cuts, in case they are needed in the future  
  new_svm_constraints <- matrix(0, nrow = 0, ncol = NVars)
  new_svm_rhs <- c()
  
  J <- 0
  Gamma <- 0 #approximator of J
  
  alpha <- NULL
  b <- 0
  eta <- NULL
  
  initial_cuts_added <- 0
  
  if(is.null(initial_info) == FALSE)
  {
    if(!is.null(initial_info$svm_cuts))
    {
      master_constraint_matrix <- rbind(master_constraint_matrix, initial_info$svm_cuts$matrix[,1:NVars])
      master_rhs <- c(master_rhs, initial_info$svm_cuts$rhs)
      master_senses <- c(master_senses, rep("G", length(initial_info$svm_cuts$rhs)))
      initial_cuts_added <- length(initial_info$svm_cuts$rhs)
      print(sprintf("Added %d initial SVM cuts.", initial_cuts_added))
    }
    
    if(!is.null(initial_info$solution))
    {
      eta <- initial_info$solution$eta
      if(G > 0)
        eta <- sparsify_eta(eta, G)
      
      b <- initial_info$solution$b
      alpha <- initial_info$solution$alpha
      J <- initial_info$solution$J
      
      incumbent_solution <- list(UB = J, eta = eta, alpha = alpha, b = b)
    }
  }
  
  objectives <- c()
  
  info_cols <- c("obj", "UB", "LB", "Gap", "Gamma", "J", "Time(M)", "Time(SVM)")
  iteration_info <- matrix(0, 0, length(info_cols))
  colnames(iteration_info) <- info_cols
  
  total_time <- 0
  
  while(terminate == FALSE)
  {
    #SOLVE MP and obtain Gamma and eta
    if(iteration > 1 | is.null(eta) | initial_cuts_added > 0)
    {
      master_result <- solve_master_problem(master_obj_coef, master_constraint_matrix, master_rhs, master_senses, 
                                                   master_lb, master_ub, master_vtype, P, contains_g = FALSE)
      Gamma <- master_result$Gamma
      eta <- master_result$eta
      master_time <- master_result$Time
    }else
    {
      Gamma <- 0
      master_time <- 0
    }
    
    eta_mp <- eta
    if(G > 0)
      eta <- sparsify_eta(eta, G)
    
    Keta <- calculate_Keta(Km, eta)
    svm_result <- solve_classification_svm(Keta, y, parameters$C, parameters$epsilon)
    
    alpha <- svm_result$alpha #regularized alpha (may be numerically infeasible)
    b <- svm_result$b
    
    alpha_cut <- svm_result$alpha_original #numerically feasible
    J <- svm_result$objective_original
    
    #model$alpha_original is already multiplied by y
    alpha_sum <- t(alpha_cut) %*% y
    kappa <- calucalte_kappa(alpha_cut, Km, P)
    
    #add constraint Gamma + Eta*kappa >= alpha_sum
    constraint_alpha <- numeric(NVars)
    constraint_alpha[1] <- 1
    constraint_alpha[2:(P+1)] <- kappa
    
    constraint_satisfied <- sum(constraint_alpha*c(Gamma, eta_mp)) >= alpha_sum
    
    if(constraint_satisfied == FALSE)
    {
      master_constraint_matrix <- rbind(master_constraint_matrix, constraint_alpha)
      master_rhs <- c(master_rhs, alpha_sum)
      master_senses <- c(master_senses, "G")  
      
      new_svm_constraints <- rbind(new_svm_constraints, constraint_alpha)
      new_svm_rhs <- c(new_svm_rhs, alpha_sum)
    }
    
    svm_time <- svm_result$Time
    
    total_time <- total_time + master_time + svm_time
    
    obj <- J
    if(obj < incumbent_solution$UB)
    {
      incumbent_solution = list(UB = obj, eta = eta, alpha = alpha, b = b)
    }
    
    LB <- Gamma
    
    Gap <- (incumbent_solution$UB-LB)/max(abs(LB),parameters$epsilon)
    
    objectives <- c(objectives, obj)
    
    info <- c(obj, incumbent_solution$UB, LB, Gap, Gamma, J, master_time, svm_time)
    iteration_info <- rbind(iteration_info, info)
    
    if(iteration %% 5 == 1)
      print(sprintf("%d %f %f %f %f%% %f %f %f %f",iteration, obj, incumbent_solution$UB, LB, Gap*100, Gamma, J, master_time, svm_time))
    
    if(constraint_satisfied | Gap <= parameters$optimality_gap | iteration >= parameters$iteration_count)
    {
      terminate <- TRUE
    }
    
    if(terminate & iteration %% 5 != 1)
      print(sprintf("%d %f %f %f %f%% %f %f %f %f",iteration, obj, incumbent_solution$UB, LB, Gap*100, Gamma, J, master_time, svm_time))
    
    iteration <- iteration + 1
  }
  
  reta <- polish_eta(incumbent_solution$eta, parameters$epsilon)
  
  state <- list(alpha = incumbent_solution$alpha, b = incumbent_solution$b, eta = reta, UB = incumbent_solution$UB,
                objectives = objectives, iteration_info = iteration_info, parameters = parameters, total_time = total_time)
  
  new_cuts_info = list(SVM_cuts = list(matrix = new_svm_constraints, rhs = new_svm_rhs))
  
  output <- list(state = state, new_cuts_info = new_cuts_info)
  
  return(output)
}

spmkl_step1_phase13_iterative <- function(Km, y, parameters) {
  G <- 0
  if(parameters$sparsify == TRUE)
    G <- parameters$G
  
  P <- dim(Km)[3]
  eta <- rep(1 / P, P)
  Keta <- calculate_Keta(Km, eta)
  
  start_time = Sys.time()
  
  model <- solve_classification_svm(Keta, y, parameters$C, parameters$epsilon)
  mkl_objectives <- model$objective
  
  incumbent_sol <- list(UB = model$objective, eta = eta, b = model$b, alpha = model$alpha)
  
  print(sprintf("******** %d: %f",length(mkl_objectives), model$objective))
  
  incumbent_sparisfied_sol <- list(UB = Inf)
  
  NVars <- 1 + P
  new_svm_constraints <- matrix(0, nrow = 0, ncol = NVars)
  new_svm_rhs <- c()
  new_svm_type <- c()
  
  bmkl_objectives <- c()
  
  while(1) {
    for (m in 1:P) {
      eta[m] <- eta[m] * sqrt(t(model$alpha) %*% Km[,,m] %*% model$alpha)
    }
    eta <- eta / sum(eta)
    eta[eta < parameters$epsilon] <- 0
    eta <- eta / sum(eta)
    Keta <- calculate_Keta(Km, eta)
    
    model <- solve_classification_svm(Keta, y, parameters$C, parameters$epsilon)
    
    if(model$objective < incumbent_sol$UB)
      incumbent_sol <- list(UB = model$objective, eta = eta, b = model$b, alpha = model$alpha)
    
    alpha_sum <- t(model$alpha_original) %*% y
    kappa <- calucalte_kappa(model$alpha_original, Km, P)
    constraint_alpha <- numeric(NVars)
    constraint_alpha[1] <- 1
    constraint_alpha[2:(P+1)] <- kappa
    new_svm_constraints <- rbind(new_svm_constraints, constraint_alpha)
    new_svm_rhs <- c(new_svm_rhs, alpha_sum)
    new_svm_type <- c(new_svm_type, "C")
    
    mkl_objectives <- c(mkl_objectives, model$objective)
    
    iteration  = length(mkl_objectives)
    
    if(iteration %% 5 == 0)
    {
      print(sprintf("******** %d: %f %f",iteration, model$objective, incumbent_sol$UB))
    }
    
    if (iteration == parameters$iteration_count) {
      parameters$stop_condition = "iteration"
      break
    }
  }
  
  if(G > 0)
  {
    eta_G <- sparsify_eta(incumbent_sol$eta, G)
    if(max(abs(eta_G-incumbent_sol$eta)) <= parameters$epsilon)
    {
      alpha = incumbent_sol$alpha
      b = incumbent_sol$b
      J = incumbent_sol$UB
    }
    else
    {
      Keta <- calculate_Keta(Km, eta_G)
      svm_result <- solve_classification_svm(Keta, y, parameters$C, parameters$epsilon)
      
      alpha <- svm_result$alpha
      b <- svm_result$b
      
      alpha_cut <- svm_result$alpha_original
      J <- svm_result$objective_original
      
      alpha_sum <- t(alpha_cut) %*% y
      kappa <- calucalte_kappa(alpha_cut, Km, P)
      constraint_alpha <- numeric(NVars)
      constraint_alpha[1] <- 1
      constraint_alpha[2:(P+1)] <- kappa
      new_svm_constraints <- rbind(new_svm_constraints, constraint_alpha)
      new_svm_rhs <- c(new_svm_rhs, alpha_sum)
      new_svm_type <- c(new_svm_type, "B")
    }
    
    bmkl_objectives <- c(bmkl_objectives, J)
    
    if(J < incumbent_sparisfied_sol$UB)
    {
      incumbent_sparisfied_sol <- list(UB = J, eta = eta_G, b = b, alpha = alpha)
    }
    
    print(sprintf("******** %d: %f %f %f %f",iteration, model$objective, J, incumbent_sol$UB, incumbent_sparisfied_sol$UB))
  }
  
  if(G > 0)
  {
    state <- list(alpha = incumbent_sparisfied_sol$alpha, b = incumbent_sparisfied_sol$b, eta = incumbent_sparisfied_sol$eta, UB = incumbent_sparisfied_sol$UB, 
                  objectives = bmkl_objectives, parameters = parameters, mkl_sol = incumbent_sol, mkl_objectives = mkl_objectives)
  }
  else
  {
    state <- list(alpha = incumbent_sol$alpha, b = incumbent_sol$b, eta = incumbent_sol$eta, UB = incumbent_sol$UB, mkl_objectives = mkl_objectives, parameters = parameters)
  }
  new_cuts_info = list(SVM_cuts = list(matrix = new_svm_constraints, rhs = new_svm_rhs, type = new_svm_type))
  
  output <- list(state = state, new_cuts_info = new_cuts_info)
  
  return(output)
}

spmkl_step1_three_phase <- function(Km, y, parameters, initial_info = NULL) {
  
  P <- dim(Km)[3]
  
  states <- list()
  final_state <- NULL
  phase_parameters <- parameters
  
  phase_2_initial_info <- initial_info
  if(!is.null(initial_info))
    initial_solution <- initial_info$solution
  
  if(parameters$phase_iterations[1] > 0)
  {
    phase_parameters$iteration_count <- parameters$phase_iterations[1]
    phase_parameters$sparsify <- TRUE
   
    if(parameters$phase1_method == "Cutting-plane")
    {
      output <- spmkl_step1_phase13_cutting_plane(Km, y, phase_parameters)
    }
    else
    {
      output <- spmkl_step1_phase13_iterative(Km, y, phase_parameters)
    }
    
    initial_solution <- list(eta = output$state$eta, g = 1*(output$state$eta>0),
                             b = output$state$b, alpha = output$state$alpha, J = output$state$UB)
    
    states[["phase1"]] <- output$state
    final_state <- output$state
    
    phase_2_initial_info <- list()
    
    phase_2_initial_info$solution <- initial_solution
    phase_2_initial_info$svm_cuts <- output$new_cuts_info$SVM_cuts
  }
  
  if(parameters$phase_iterations[2] > 0)
  {
    phase_parameters$iteration_count <- parameters$phase_iterations[2]
    phase_parameters$sparsify <- FALSE
    
    output <- spmkl_step1_phase2_cutting_plane(Km, y, phase_parameters, initial_info = phase_2_initial_info)
    final_state <- output$state
    initial_solution <- list(eta = output$state$eta, g = output$state$g,
                             b = output$state$b, alpha = output$state$alpha, J = output$state$UB)
    
    states[["phase2"]] <- output$state
  }
  
  if(parameters$phase_iterations[3] > 0)
  {
    phase_parameters$iteration_count <- parameters$phase_iterations[3]
    phase_parameters$sparsify <- FALSE
    
    selected_kernels <- which(initial_solution$g > 0)
    phase3_Km <- Km[,,selected_kernels]
    
    new_P <- length(selected_kernels)
    phase3_state <- NULL
    if(new_P > 1)
    {
      if(parameters$phase3_method == "Cutting-plane")
      {
        phase_3_initial_info <- list(solution = initial_solution)
        phase_3_initial_info$solution$eta <- phase_3_initial_info$solution$eta[selected_kernels]
        output <- spmkl_step1_phase13_cutting_plane(phase3_Km, y, phase_parameters, initial_info = phase_3_initial_info)
      }
      else
      {
        output <- spmkl_step1_phase13_iterative(phase3_Km, y, phase_parameters)
      }
      
      phase3_state <- output$state
    }
    else
    {
      phase3_state <- list(eta = 1, g = initial_solution$g, b = initial_solution$b, alpha = initial_solution$alpha, 
                           UB = initial_solution$J, total_time = 0, parameters = phase_parameters)
    }
    
    new_eta <- numeric(P)
    new_eta[selected_kernels] <- phase3_state$eta
    phase3_state$eta <- new_eta
    
    final_state <- phase3_state
    states[["phase3"]] <- phase3_state
  }
  
  output <- list(state = final_state, states = states, parameters = parameters)
  
  return(output)
}

spmkl_step2_iterative <- function(Km, y, parameters, kernel_significance = NULL) {
  
  get_eta_error <- function(et)
  {
    return((sum(et)-1)^2)
  }
  
  compute_eta <- function(f, t, error_tol = 1e-10, max_reps = 100, print_results = FALSE)
  {
    et <- f/sum(f)
    best_eta <- et
    min_error <- Inf
    
    if(is.null(t) == FALSE)
    {
      for(r in 1:max_reps)
      {
        lambda_j <- t+(f/et)^2
        lambda <- min(lambda_j[which(lambda_j > max(t))])
        et <- f/sqrt(lambda-t)
        error <- get_eta_error(et)
        
        et <- et / sum(et)
        
        if(error < min_error)
        {
          min_error <- error
          best_eta <- et
        }
        
        if(error < error_tol)
        {
          break;
        }
      }
    }
    
    return(best_eta)
  }
  
  theta <- NULL
  if(is.null(kernel_significance) == FALSE)
    theta <- kernel_significance * parameters$penalty_weight
  
  P <- dim(Km)[3]
  eta <- rep(1 / P, P)
  Keta <- calculate_Keta(Km, eta)
  
  start_time = Sys.time()
  
  model <- solve_classification_svm(Keta, y, parameters$C, parameters$epsilon)
  
  objective <- model$objective_original
  if(is.null(kernel_significance) == FALSE)
    objective <- parameters$penalty_weight + model$objective_original -sum(eta*theta)
  
  mkl_objectives <- objective
  
  incumbent_sol <- list(UB = objective, eta = eta, b = model$b, alpha = model$alpha)
  
  print(sprintf("******** %d: %f",length(mkl_objectives), objective))
  
  iteration <- 1
  
  while(1) {
    
    report_info <- iteration %% 5 == 0
    
    f_tilde <- numeric(P)
    for (m in 1:P) {
      f_tilde[m] <- eta[m] * sqrt(t(model$alpha) %*% Km[,,m] %*% model$alpha)
    }
    eta <- compute_eta(f_tilde, theta, print_results = report_info)
    eta[eta < parameters$epsilon] <- 0
    eta <- eta / sum(eta)
    Keta <- calculate_Keta(Km, eta)
    
    model <- solve_classification_svm(Keta, y, parameters$C, parameters$epsilon)
    
    objective <- model$objective_original
    if(is.null(kernel_significance) == FALSE)
      objective <- parameters$penalty_weight + model$objective_original -sum(eta*theta)
    
    if(objective < incumbent_sol$UB)
      incumbent_sol <- list(UB = objective, eta = eta, b = model$b, alpha = model$alpha)
    
    mkl_objectives <- c(mkl_objectives, objective)
    
    if(report_info)
    {
      print(sprintf("******** %d: %f %f",iteration, objective, incumbent_sol$UB))
    }
    
    if (iteration == parameters$iteration_count) {
      parameters$stop_condition = "iteration"
      
      if(!report_info)
      {
        print(sprintf("******** %d: %f %f",iteration, objective, incumbent_sol$UB))
      }
      
      break
    }
    
    iteration <- iteration + 1
  }
  
  state <- list(alpha = incumbent_sol$alpha, b = incumbent_sol$b, eta = incumbent_sol$eta, 
                UB = incumbent_sol$UB, mkl_objectives = mkl_objectives, parameters = parameters)
  output <- list(state = state)
  
  return(output)
}
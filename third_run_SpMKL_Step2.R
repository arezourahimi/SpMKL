library(AUC)
code_path <- "./code"
pathway <- 'hallmark'

source(sprintf("%s/spmkl_train.R", code_path))
source(sprintf("%s/spmkl_test.R", code_path))
source(sprintf("%s/solve_classification_models_cplex.R", code_path))
source(sprintf("%s/classification_helper.R", code_path))

cohorts <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-ESCA", "TCGA-HNSC", "TCGA-KICH",
             "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC",
             "TCGA-PAAD", "TCGA-READ", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA")

TN <- length(cohorts)

data_path <- "../data"
result_path <- "./results_SpMKL_step2"
if (dir.exists(result_path) == FALSE) {
  dir.create(result_path)
}

kernel_significance_file <- sprintf("%s/kernel_significance_matrix.csv",data_path)
if(file.exists(kernel_significance_file) == FALSE)
  stop("kernel signficance matrix not found!")
kernel_significance <- read.csv(kernel_significance_file, row.names = 1)
colnames(kernel_significance) <- cohorts

pathways <- read_pathways(pathway)
gene_names <- sort(unique(unlist(sapply(1:length(pathways), FUN = function(x) {pathways[[x]]$symbols}))))

X <- vector("list", TN)
y <- vector("list", TN)
negative_indices <- vector("list", TN)
positive_indices <- vector("list", TN)
for (t in 1:TN) {
  load(sprintf("%s/%s.RData", data_path, cohorts[t]))
  
  common_patients <- intersect(rownames(TCGA$clinical)[which(is.na(TCGA$clinical$pathologic_stage) == FALSE)], rownames(TCGA$mrna))
  
  X[[t]] <- log2(TCGA$mrna[common_patients,] + 1)
  y[[t]] <- rep(NA, length(common_patients))
  
  if(stage == "1vs234")
  {
    y[[t]][TCGA$clinical[common_patients, "pathologic_stage"] %in% c("Stage I",  "Stage IA",  "Stage IB",  "Stage IC")] <- +1
    y[[t]][TCGA$clinical[common_patients, "pathologic_stage"] %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC",
                                                                     "Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC",
                                                                     "Stage IV",  "Stage IVA",  "Stage IVB",  "Stage IVC")] <- -1
  }else
  {
    y[[t]][TCGA$clinical[common_patients, "pathologic_stage"] %in% c("Stage I",  "Stage IA",  "Stage IB",  "Stage IC",
                                                                     "Stage II", "Stage IIA", "Stage IIB", "Stage IIC")] <- +1
    y[[t]][TCGA$clinical[common_patients, "pathologic_stage"] %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC",
                                                                     "Stage IV",  "Stage IVA",  "Stage IVB",  "Stage IVC")] <- -1
  }
  
  valid_patients <- which(is.na(y[[t]]) == FALSE)
  valid_features <- as.numeric(which(apply(X[[t]][valid_patients,], 2, sd) != 0))
  X[[t]] <- X[[t]][valid_patients, valid_features]
  y[[t]] <- y[[t]][valid_patients]
  
  negative_indices[[t]] <- which(y[[t]] == -1)
  positive_indices[[t]] <- which(y[[t]] == +1)
  
  X[[t]] <- X[[t]][, which(colnames(X[[t]]) %in% gene_names)]
}

P <- length(pathways)

base_seed <- 1505

replication_count <- 100
replications <- 1:replication_count
cv_tuples <- get_cross_validation_tuples_step2()

for (t in 1:TN) {
  cohort <- cohorts[t]
  
  if (dir.exists(sprintf("%s/%s", result_path, cohort)) == FALSE) {
    dir.create(sprintf("%s/%s", result_path, cohort)) 
  }
  
  for(replication in replications)
  {
    state_file <- sprintf("%s/%s/spmkl_step2_rep_%d_state.RData", result_path, cohort, replication)
    
    if (file.exists(state_file) == FALSE){
      
      epsilon <- 1e-5
      fold_count <- 4
      train_ratio <- 0.8
      iteration_count <- 200
      
      set.seed(base_seed * replication)
      train_negative_indices <- sample(negative_indices[[t]], ceiling(train_ratio * length(negative_indices[[t]])))
      train_positive_indices <- sample(positive_indices[[t]], ceiling(train_ratio * length(positive_indices[[t]])))
      
      auroc_tuples <- matrix(0, nrow = nrow(tuples), ncol = 2+fold_count)
      colnames(auroc_tuples) <- c("C", "Penalty", paste("Auroc",1:fold_count))
      
      negative_allocation <- sample(rep(1:fold_count, ceiling(length(train_negative_indices) / fold_count)), length(train_negative_indices))
      positive_allocation <- sample(rep(1:fold_count, ceiling(length(train_positive_indices) / fold_count)), length(train_positive_indices))
      
      cv_start_time <- Sys.time()
      
      for (fold in 1:fold_count) {
        train_indices <- c(train_negative_indices[which(negative_allocation != fold)], train_positive_indices[which(positive_allocation != fold)])
        test_indices <- c(train_negative_indices[which(negative_allocation == fold)], train_positive_indices[which(positive_allocation == fold)])
        
        X_train <- X[[t]][train_indices,]
        X_test <- X[[t]][test_indices,]
        X_train <- scale(X_train)
        X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
        
        N_train <- nrow(X_train)
        N_test <- nrow(X_test)
        y_train <- y[[t]][train_indices]
        y_test <- y[[t]][test_indices]
        
        K_train <- array(0, dim = c(N_train, N_train, P))
        K_test <- array(0, dim = c(N_test, N_train, P))
        for (m in 1:P) {
          feature_indices <- which(colnames(X_train) %in% pathways[[m]]$symbols)
          D_train <- pdist(X_train[, feature_indices], X_train[, feature_indices])
          D_test <- pdist(X_test[, feature_indices], X_train[, feature_indices])
          sigma <- mean(D_train)
          K_train[,,m] <- exp(-D_train^2 / (2 * sigma^2))
          K_test[,,m] <- exp(-D_test^2 / (2 * sigma^2))
        }
        
        for (tpl in 1:nrow(cv_tuples)) {
          C_tpl <- tuples[tpl, "C"]
          penalty_tpl <- tuples[tpl, "Penalty"]
          
          print(sprintf("%s Replication %d --> fold = %d, C = %g, Penalty = %g", cohort, replication, fold, C_tpl, penalty_tpl))
          parameters <- list()
          parameters$C <- C_tpl
          parameters$penalty_weight <- penalty_tpl
          parameters$epsilon <- epsilon
          parameters$iteration_count <- iteration_count
          
          output <- spmkl_step2_iterative(K_train, y_train, parameters, kernel_significance)
          state <- output$state
          prediction <- spmkl_test(K_test, state)
          AUROC <- auc(roc(prediction$f, as.factor(y_test)))
          
          auroc_tuples[tpl, c("C", "Penalty", paste("Auroc", fold))] <- c(C_tpl, penalty_tpl, AUROC)
        }
      }
      
      cv_total_time <- as.numeric(Sys.time()-cv_start_time, unit = "secs")
      
      test_start_time <- Sys.time()
      
      average_aurocs <- rowMeans(auroc_tuples[,2+(1:fold_count)])
      tuple_star <- which(average_aurocs == max(average_aurocs))[1]
      C_star <- auroc_tuples[tuple_star, "C"]
      Penalty_star <- auroc_tuples[tuple_star, "Penalty"]
      
      train_indices <- c(train_negative_indices, train_positive_indices)
      test_indices <- setdiff(1:length(y[[t]]), train_indices)
      
      X_train <- X[[t]][train_indices,]
      X_test <- X[[t]][test_indices,]
      X_train <- scale(X_train)
      X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
      
      N_train <- nrow(X_train)
      N_test <- nrow(X_test)
      
      K_train <- array(0, dim = c(N_train, N_train, P))
      K_test <- array(0, dim = c(N_test, N_train, P))
      for (m in 1:P) {
        feature_indices <- which(colnames(X_train) %in% pathways[[m]]$symbols)
        D_train <- pdist(X_train[, feature_indices], X_train[, feature_indices])
        D_test <- pdist(X_test[, feature_indices], X_train[, feature_indices])
        sigma <- mean(D_train)
        K_train[,,m] <- exp(-D_train^2 / (2 * sigma^2))
        K_test[,,m] <- exp(-D_test^2 / (2 * sigma^2))
      }
      
      y_train <- y[[t]][train_indices]
      y_test <- y[[t]][test_indices]
      
      parameters <- list()
      parameters$C <- C_star
      parameters$penalty_weight <- Penalty_star
      parameters$epsilon <- epsilon
      parameters$iteration_count <- iteration_count
      
      output <- spmkl_step2_iterative(K_train, y_train, parameters, kernel_significance)
      state <- output$state
      state$auroc_cv_tuples <- auroc_cv_tuples
      prediction <- spmkl_test(K_test, state)
      state$AUROC <- auc(roc(prediction$f, as.factor(y_test)))
      state$test_time <- as.numeric(Sys.time()-test_start_time, unit = "secs")
      state$cv_time <- cv_total_time
      
      save("state", file = state_file)
    } 
  }
}

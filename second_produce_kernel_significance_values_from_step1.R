cohorts <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-ESCA", "TCGA-HNSC", "TCGA-KICH",
             "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC",
             "TCGA-PAAD", "TCGA-READ", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA")
TN <- length(cohorts)

replication_count <- 100
replications <- 1:replication_count

data_path <- "../data"
result_path <- "./results_SpMKL_step1"
kernel_significance_file <- sprintf("%s/kernel_significance_matrix.csv", data_path)

pathway <- 'hallmark'
pathways <- read_pathways(pathway)
pathway_names <- sapply(pathways, function(p) {p$name})
if(startsWith(pathway_names[1], "HALLMARK"))
  pathway_names <- substr(pathway_names, nchar("HALLMARK") + 2, nchar(pathway_names))

P <- length(pathways)

kernel_significance_matrix <- matrix(0, nrow = P, ncol = TN, dimnames = list(pathway_names, cohorts))

for (t in 1:TN) {
  for(replication in replications)
  {
    state_file <- sprintf("%s/%s/spmkl_step1_rep_%d_state.RData", result_path, cohorts[t], replication)
    load(state_file)
    kernel_significance_matrix[, t] <- state$g
  }
}

kernel_significance_matrix <- kernel_significance_matrix / replication_count

write.csv(kernel_significance_matrix, file = kernel_significance_file)
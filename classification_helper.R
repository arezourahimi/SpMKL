pdist <- function(X1, X2) {
  if (identical(X1, X2) == TRUE) {
    D <- as.matrix(dist(X1))
  }
  else {
    D <- as.matrix(dist(rbind(X1, X2)))
    D <- D[1:nrow(X1), (nrow(X1) + 1):(nrow(X1) + nrow(X2))]
  }
  return(D)
}

read_pathways <- function(name) {
  symbols_lines <- read.table(sprintf("msigdb/%s.gmt", name), header = FALSE, sep = ",", stringsAsFactor = FALSE)
  pathways <- vector("list", nrow(symbols_lines))
  for (line in 1:nrow(symbols_lines)) {
    symbols_entries <- strsplit(symbols_lines[line, 1], "\t")
    pathways[[line]]$name <- symbols_entries[[1]][1]
    pathways[[line]]$link <- symbols_entries[[1]][2]
    pathways[[line]]$symbols <- sort(symbols_entries[[1]][-2:-1])
  }
  return(pathways)
}

calculate_Keta <- function(Km, eta) {
  P <- dim(Km)[3]
  
  Keta <- eta[1] * Km[,,1]
  for (m in 2:P) {
    Keta <- Keta + eta[m] * Km[,,m]
  }
  
  return(Keta)
}

#Calculates the coefficients of kernel weights in the objective function of dual SVM problem for a given alpha
calucalte_kappa <- function(alpha, Km, P)
{
  kappa <- numeric(P)
  
  for(m in 1:P)
  {
    kappa[m] <- 0.5*t(alpha) %*% Km[,,m] %*% alpha
  }
  return(kappa)
}

#Removes the numerically insignificat kernel weights
polish_eta <- function(eta, epsilon)
{
  et <- eta / sum(eta)
  et[et < epsilon] <- 0
  et <- et / sum(et)
  return(et)
}

sparsify_eta <- function(eta, G)
{
  ord <- order(eta, decreasing = T)[1:G]
  et <- numeric(length(eta))
  et[ord] <- eta[ord]
  et <- et / sum(et)
  return(et)
}

#Yields different combinations of C and G
get_cross_validation_tuples_step1 <- function()
{
  C_set <- 10^(-4:5)
  G_set <- c(5, 10, 15, 20)
  
  tuples <- as.matrix(expand.grid(C_set, G_set))
  colnames(tuples) <- c("C", "G")
  
  return(tuples)
}

#Yields different combinations of C and rho
get_cross_validation_tuples_step2 <- function()
{
  C_set <- C
  rho_set <- 10^(-1:3)
  
  tuples <- as.matrix(expand.grid(C_set, rho_set))
  colnames(tuples) <- c("C", "Penalty")
  
  return(tuples)
}



# X : Matrix of molecular markers coded as 0 (A1A1), 1 (A1A2) and 2 (A2A2); or matrix of envirotype data (mean scaled and centered)
# is.center = FALSE means that X is not mean centered. If X is already mean centered, use is.center = TRUE
# output: Linear covariance matrix

GB_Kernel <-function(X,is.center=FALSE){
  if(isFALSE(is.center)) X = scale(x = X,center = T,scale = F)
  XXl <- X %*% t(X)
  K_G <- XXl/(sum(diag(XXl))/nrow(X)) + diag(1e-6, nrow(XXl))
  return(K_G)
}

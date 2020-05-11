# M : Matrix of molecular markers coded as 0 (A1A1), 1 (A1A2) and 2 (A2A2)
# output: Ga = Genomic Relationship Matrix (GRM) for additive effects
#       : Gd = GRM for Dominance Effects
#       : M = returns the original M matrix
#       : S = returns de molecular matrix re-codes for dominance deviaitons

GB_kernel <-function(M){
  N <- nrow(M)
  m <- ncol(M)
  p <- colMeans(M)/2
  
  WWG <- function(M, p){
    w <- scale(x = M, center = T, scale = F)
    
    S <- ((M==2)*1) * - rep(2*(1-p)^2, each=N) + ((M==1)*1) * rep(2*p*(1-p), each=N) + ((M==0)*1) * (-rep(2*p^2, each=N))
    
    WWl <- w %*% t(w)
    Ga <- WWl/(sum(diag(WWl))/N) + diag(1e-6, nrow(WWl))
    
    SSl <- S %*% t(S)
    Gd <- SSl/(sum(diag(SSl))/N)
    
    return(list(Ga=Ga,Gd=Gd,M=M,S=S))
  }
  
  return(WWG(M=M,p=p))
}

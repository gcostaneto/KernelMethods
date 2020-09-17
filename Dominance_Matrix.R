Dominance <-function(M){
  N <- nrow(M)
  m <- ncol(M)
  p <- colMeans(M)/2
  
  WWG <- function(M, p){
    w <- scale(x = M, center = T, scale = F)
    
    S <- ((M==2)*1) * - rep(2*(1-p)^2, each=N) + ((M==1)*1) * rep(2*p*(1-p), each=N) + ((M==0)*1) * (-rep(2*p^2, each=N))
    
    return(S)
  }
  
  return(WWG(M=M,p=p))
}

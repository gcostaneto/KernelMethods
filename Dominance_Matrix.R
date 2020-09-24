# Dominance effects matrix (S)
# M is a molecular marker matrix (0,1,2)

S.matrix <-function(M){
  # from SnpReady
  N <- nrow(M)
  m <- ncol(M)
  p <- colMeans(M)/2
  w <- scale(x = M, center = T, scale = F)
  S <- ((M==2)*1) * - rep(2*(1-p)^2, each=N) + ((M==1)*1) * rep(2*p*(1-p), each=N) + ((M==0)*1) * (-rep(2*p^2, each=N))
 return(S) 
}

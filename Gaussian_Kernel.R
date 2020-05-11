# X = list of molecular or environmental data matrix
# h = bandwith factor (constant)
# prob = quantile used for GK (default is the 2th quantile = 0.5 = median)
# y = phenotypic data. If y = NULL, the bandwith used is h. Else, a marginal function are used to compute the better banddwith factor

GK_Kernel <- function (X, h = NULL,prob=.5,y=NULL) 
{
  
  margh.fun <- function(theta, y, D, q, nu=0.0001, Sc=0.0001, nuh=NULL, Sch=NULL, prior=NULL)
  {
    h <- theta[1]
    phi <- theta[2]
    Kh <- exp(-h*D/q)
    eigenKh <- eigen(Kh)
    nr <- length(which(eigenKh$val> 1e-10))
    Uh <- eigenKh$vec[,1:nr]
    Sh <- eigenKh$val[1:nr]
    d <- t(Uh) %*% scale(y, scale=F)
    
    Iden <- -1/2*sum(log(1+phi*Sh)) - (nu+nr-1)/2*log(Sc+sum(d^2/(1+phi*Sh)))
    if(!is.null(prior)) lprior <- dgamma(h,nuh,Sch,log=T) else lprior <- 0
    
    Iden <- -(Iden+lprior)
    
    return(Iden)
  }
  
  GK <-list()
  for(i in 1:length(X)){
    d <- as.matrix(dist(X[[i]], upper = T, diag = T))^2
    q <- quantile(x=d,prob=prob)
    if(q == 0) q <- quantile(x=d,prob=.05)
    if (is.null(h)) h <- 1
    
    if(!is.null(y)){
      sol<-optim(c(1,3),margh.fun,y=y,D=d,q=q,method="L-BFGS-B",
                 lower=c(0.05,0.05),upper=c(6,30),prior="gamma",nuh=3,Sch=1.5)
      h<-sol$par[1]
    }
    
    
    GK[[i]] <- exp(-h * d/q)
  }
  names(GK) <- names(X)
  return(GK)
}

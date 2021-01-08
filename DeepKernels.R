#'##############################################################################
#' Deep kernel functions based on Arc-cosine kernels
#' Authors: Jaime Cuevas (DK) & Germano Costa Neto (DK for reaction-norm)
#'##############################################################################
# last update: Jan 08 2020

#'------------------------------------------------------------------------
# Create DK for each Matrix (M) of inputs AK (Author: GCN) #####################
#'------------------------------------------------------------------------

# same as GC1.fun, but for multiples elements into a M list (ex: M = list(Additive,Dominance)) -------------------------------------
get_GC1 <- function(M){
  AK1 <- list()
  for(i in 1:length(M)) AK1[[i]] <- GC1.fun(X = M[[i]])
  length(AK1)
  names(AK1) = names(M)
  return(AK1)
  
}

#'------------------------------------------------------------------------
# Optimzation of DK (integrating marginal functions) #####################
#'------------------------------------------------------------------------
#' K  : kernel output from get_kernel function of EnvRtype
#' y  : phenotypic records (with NAs)
#' tr : training set identification
#' nl : predeterminated maximum number of hidden layers
#' package: if you want to run in another package, use package = 'other'. If you wanto run in BGGE or EnvRtype, package = 'BGGE' by default

opt_AK <- function(K,y, tr, nl=40,package = 'BGGE')
  {
  id <- names(K)
  .K_post <-list()
  if(package == 'BGGE'){for(j in 1:length(K)) .K_post[[j]] <- K[[j]]$Kernel} 
  if(!package == 'BGGE'){for(j in 1:length(K)) .K_post[[j]] <- K[[j]]} 
  opt_K  <- list()
  for(i in 1:length(K))
    {
    l          <- marg.AK(y=y,GC=.K_post[[i]][tr,tr], nl=nl)
    cat(paste0(Sys.time(),'  Deep Kernel for: ',id[i],' effect with ',l, ' layers \n'))
    opt_K[[i]] <- Kernel.function(GC=.K_post[[i]],nl=l)
  }
  if(package == 'BGGE'){for(j in 1:length(K)) K[[j]]$Kernel <- opt_K[[i]]} 
  if(!package == 'BGGE'){for(j in 1:length(K)) K[[j]]<- opt_K[[i]]} 
  
  return(K)
}


#'------------------------------------------------------------------------
# Marginal function for AK (Author: JC) ########################################
#'------------------------------------------------------------------------
marg.AK <- function(y,GC,nl){
  lden.fun<-function(theta,nr,Uh,Sh,d){
    phi<-theta[1]
    lden  <- -1/2*sum(log((1+phi*Sh)))-(nr-1)/2*log(sum(d^2/((1+phi*Sh))))
    lden <- -(lden)
    return(lden)
  }
  vero<-function(y,GC) {          
    Kh <- GC
    eigenKh <- eigen(Kh)
    nr<- length(which(eigenKh$val>1e-10))
    Uh <- eigenKh$vec[,1:nr]
    Sh <- eigenKh$val[1:nr]
    d <- t(Uh)%*%scale(y,scale=F)
    sol <-optim(c(1),lden.fun,nr=nr,Uh=Uh,Sh=Sh,d=d,method="L-BFGS-B",lower=c(0.0005),upper=c(200))
    phi<-sol$par[1]
    log.vero<--1/2*sum(log((1+phi*Sh)))-(nr-1)/2*log(sum(d^2/((1+phi*Sh))))
    return(log.vero)
  }
  l<-1
  GC2<-GC
  vero1<-vero(y=y,GC=GC2)
  m<-0
  while( m==0 && (l<nl)){
    l<-l+1
    GC<-Kernel.function(GC=GC2,nl=1)
    GC2<-GC
    vero2<-vero(y=y,GC=GC2)
    if(vero2<vero1) m=1
    vero1<-vero2
  }
  return(l-1)
}


#'------------------------------------------------------------------------
# Base AK kernel  (Author: JC) #################################################
#'------------------------------------------------------------------------
GC1.fun<-function(X){
  n<-nrow(X)
  cosalfa<-cor(t(X))
  angulo<-acos(cosalfa)
  mag<-sqrt(apply(X,1,function(x) crossprod(x)))
  sxy<-tcrossprod(mag)
  GC1<-(1/pi)*sxy*(sin(angulo)+(pi*matrix(1,n,n)-angulo)*cosalfa)
  GC1<-GC1/median(GC1)
  colnames(GC1)<-rownames(X)
  rownames(GC1)<-rownames(X)
  
  return(GC1)
}  


Kernel.function<-function(GC,nl){
  n<-nrow(GC)
  GC1<-GC
  
  for ( l in 1:nl){
    
    Aux<-tcrossprod(diag(GC))
    cosalfa<-GC*(Aux^(-1/2))
    cosa<-as.vector(cosalfa)
    cosa[which(cosalfa>1)]<-1
    
    angulo<-acos(cosa)
    angulo<-matrix(angulo,n,n)
    
    GC<-(1/pi)*(Aux^(1/2))*(sin(angulo)+(pi*matrix(1,n,n)-angulo)*cos(angulo))
    
  }
  
  GC<-GC/median(GC)
  
  rownames(GC)<-rownames(GC1)
  colnames(GC)<-colnames(GC1)
  return(GC)
}

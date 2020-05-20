#'##############################################################################
#' Deep kernel functions based on Arc-cosine kernels
#' Authors: Jaime Cuevas & Germano Costa Neto
#'##############################################################################

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

# Optimzation of DK (integrating marginal functions) -------------------------------------
# model: type of model output (same as get_kernels() function)
# MM (main effect), MDs (MM+GE), RNMM, RNMDs, see get_kernels function from Envrtype package

opt_AK <- function(K_G,K_E=NULL,nl=40,Y,tr,model='MM'){
  
  y <- Y$value[tr]
  
  #' Creating basic kernels -----------
  if(model == 'MDs' | model == 'EMDs'){
    GC <- get_kernel(K_G = K_G,K_E = K_E,model = 'MM',Y = Y)
  }else{
    if(!model == 'MDs' | !model == 'EMDs') GC <- get_kernel(K_G = K_G,K_E = K_E,model = model,Y = Y)
  } 
  
  
  gids <- GC[grep(names(GC),pattern = 'KG_' )]
  envs <- GC[grep(names(GC),pattern = 'KE_' )]
  
  
  #' Genetic  relatedness ----------
  GC <-list()
  for(j in 1:length(gids)) GC[[j]] <- gids[[j]]$Kernel
  
  optAKg <- list()
  
  for(i in 1:length(GC)){
    l <- marg.AK(y=y,GC=GC[[i]][tr,tr], nl=nl)
    optAKg[[i]] <- Kernel.function(GC=K_G[[i]],nl=l)
    
  }
  names(optAKg) = names(K_G)
  cat(paste0('AK for G effects ', length(optAKg),' Kernels \n'))
  
  #' Environmental relatedness ----------
  optEK <- NULL
  # step 5: optmizing E
  if(model %in% c('EMM','EMDs','RNMM','RNMDs')){
    cat(paste0('AK for E effects \n'))
    EC <-list()
    for(j in 1:length(envs)) EC[[j]] <- envs[[j]]$Kernel
    optEK <- list()
    for(i in 1:length(EC)){
      l <-marg.AK(y=y,GC=EC[[i]][tr,tr], nl=nl)
      optEK[[i]] <- Kernel.function(GC=K_E[[i]],nl=l)
      
    }
    names(optEK) = names(K_E)
  }
  
  #' using get_kernel to create multiple genomic x enviromic kernels
  
  if(!model =='MDs'){
    K <- get_kernel(K_G = optAKg,K_E=optEK,model = model,Y = Y)
  }else{
    if(model == 'MDs'){
      K <- get_kernel(K_G = optAKg,K_E=NULL,model = 'MM',Y = Y)
      nK <-names(K)
      Ze <- model.matrix(~0+env,Y)
      ZZ <- tcrossprod(Ze)
      optGE<-list()
      for(i in 1:length(optAKg)){
        optGE[[i]] <- list(Kernel=K[[i]]$Kernel*ZZ, Type='BD')
      }
      names(optGE) <- paste0('KGE_',names(optAKg),'E')
      K <- cbind(K,optGE)
      names(K) <- c(nK, names(optGE))
    }
  }
  
  
  return(K)
}

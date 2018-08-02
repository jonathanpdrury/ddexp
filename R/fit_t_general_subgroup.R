# Fit a model for which rates depends on a time-serie curve with regime specific parameters estimates.
# Initially, this was intended to fit model to a subgroup, taking a simmap from CreateGeobyClassObject_mvMORPH
# along with a dataset, then internally trimming the VCV matrix
# however, for an unknown reason, this crashes mvLL() in some part of the parameter space visited during optimization
# Current methods for fitting subgroup are shown in 20180802_testing.R, and should be wrapped into another script


fit_t_general_subgroup <- function(full.simmap.tree, data, fun, error=NULL, beta=NULL, sigma=NULL, model=c("exponential","linear"), method=c("L-BFGS-B","BB"), upper=Inf, lower=-20, control=list(maxit=20000), diagnostic=TRUE, echo=TRUE, constraint=TRUE) {
  
  tree<-full.simmap.tree
  
  require(mvMORPH)
  if(!inherits(tree,"simmap")==TRUE) stop("For now only simmap-like mapped trees are allowed.","\n")
  
  # Parameters
  if(is.null(names(data))) {stop("data should have names matching names in tree")}
   
  data<-as.matrix(data)
  method=method[1]
  #rownames(data)<-names(data)
  model=model[1]
  # Compute node time from the root to the tips
  times<-max(nodeHeights(tree))-nodeHeights(tree)[match(1:tree$Nnode+length(tree$tip),tree$edge[,1]),1]
  names(times)<-1:tree$Nnode+length(tree$tip)
  # Set the root to zero
  times<-max(times)-times
  # Max time
  mtot=max(nodeHeights(tree))
  
  
  # Number of species
  n=length(tree$tip.label)
  # Number of traits (for future versions)
  k=1
  # Number of maps (selective regimes)
  if(constraint){
      number_maps <- 1
  }else{
      number_maps <- ncol(tree$mapped.edge)
  }
  
  # Param likelihood contrasts (we are organizing the tree to be postorder for an efficient traversal)
  ind=reorder(tree,"postorder",index.only=TRUE)
  phy=tree
  phy$edge.length<-phy$edge.length[ind]
  phy$edge<-phy$edge[ind,]
  
  # check for simmap like format
  if(inherits(tree,"simmap")){
    phy$mapped.edge<-phy$mapped.edge[ind,]
    phy$maps<-phy$maps[ind]
  }
  #phy <- reorderSimmap(tree, order="pruningwise")
  
  # Random starting value if not provided
  if(is.null(beta)){
    beta=rep(0, number_maps)
  }
  if(is.null(sigma)){
    #sigma=rep(sqrt(var(data)/max(nodeHeights(extract.clade(tree,getMRCA(tree,rownames(data)))))), number_maps)
	sigma=rep(sum(pic(data,drop.tip(tree,tree$tip.label[which(!tree$tip.label%in%rownames(data))]))^2)/n, number_maps)
  }
  
  if(model=="linear"){
    startval=c(beta,sigma)
    nbeta=length(beta)
    nsigma=length(sigma)
  }else if(model=="exponential"){
    startval=c(beta,log(sigma))
    nbeta=length(beta)
    nsigma=length(sigma)
  }
  
  # Error estimation?
  if(!is.null(error)){
    ## Index error
    index_error<-sapply(1:n, function(x){ which(phy$edge[,2]==x)})
    startval=c(startval,0.001)
  }
  
  ##--------------Fonction-generale-DD-Env-------------------------------------------##
  
  BranchtransformMAPS<-function(phy,beta,mtot,times,fun,sigma=NULL,model=NULL,errorValue=NULL){
    #Transformations
    tips <- length(phy$tip.label)
    res <- phy
    
    if(model=="exponential"){  
      # because "times" assume that the root state is at 0 and funEnv is from the past to the present.
      # this will only work in the constraint=TRUE instance
      if(constraint){
      #f<-function(x, sigma, beta, funInd){(sigma*sum(exp(beta*fun[[funInd]](x))))/length(fun[[funInd]](x))}
      #f<-function(x, sigma, beta, funInd){(sigma*rowSums(exp(beta*fun[[funInd]](x))))/ncol(fun[[funInd]](x))}
     # f<-function(x, sigma, beta, funInd){sigma*(rowSums(exp(beta*fun[[funInd]](x)))/ncol(fun[[funInd]](x)))}
        f<-function(x, sigma, beta, funInd){sigma*exp(beta*fun[[funInd]](x))}
      } else{
      f<-function(x, sigma, beta, funInd){sigma*exp(beta*fun[[funInd]](x))}
      }
    }else if(model=="linear"){
      # Clim-lin function
      if(constraint){
      f<-function(x, sigma, beta, funInd){sigma+beta*fun[[funInd]](x)}
      } else{
      f<-function(x, sigma, beta, funInd){sigma+beta*fun[[funInd]](x)}
      }
    }
    
    # Loops over the edges
    for (i in 1:length(phy$edge.length)) {
      
      age <- times[phy$edge[i, 1] - tips] # retrieve the age at the node
      currentmap<-phy$maps[[i]]           # retrieve the corresponding maps on the edge "i"
      indlength<-length(currentmap)       # How many mapping there are?
      tempedge<-numeric(indlength)        # temporary vector for the mappings
      
      # loop pour traverser les "maps"
      for(betaval in 1:indlength){
        regimenumber <- which(colnames(phy$mapped.edge)==names(currentmap)[betaval])  # retrieve the regimes within maps
        if(constraint){
            bet <- beta
            sig <- sigma
        }else{
            bet<-beta[regimenumber]           # select the corresponding parameter for beta
            sig<-sigma[regimenumber]          # select the corresponding parameter for
        }
        sigma
        bl<-currentmap[[betaval]]         # branch length under the current map
        
        int <- try(integrate(f, lower=age, upper=(age + bl), subdivisions=500, rel.tol = .Machine$double.eps^0.05, sigma=sig, beta=bet, funInd=regimenumber), silent=TRUE)
                        if(inherits(int ,'try-error')){
                          warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
                          tempbeta <- NA_real_
                        } else {
                          tempbeta <- int$value
                        }
        tempedge[betaval] <- tempbeta 
        # on met à jour age parcequ'on va passer au maps suivant pour la lignée i
        # update "age" because we're moving to the next map for lineage i.
        age<-age+bl
      }
      # update branch length
      res$edge.length[i]<-sum(tempedge)
    }
    
    phy<-res
    
    if(!is.null(errorValue)){
      phy$edge.length[index_error]<-phy$edge.length[index_error]+errorValue^2
    }
    
    return(phy)
  }
  
  
  ##---------------------------------------------------------------------------------##
  clikCLIM <- function( param, dat, phylo, mtot, times, fun=fun, model, results=FALSE) {
    
    if(model=="exponential"){
      beta<-param[seq_len(nbeta)]
      sigma<-exp(param[nbeta+seq_len(nsigma)])
      if(!is.null(error)) errorValue <- param[nbeta+nsigma+1] else errorValue <- NULL
      phylo <- BranchtransformMAPS(phylo, beta, mtot, times, fun, sigma, model, errorValue)
      
      V<-vcv(phylo)
	  itoremove<-which(!colnames(V)%in%rownames(dat))
	  if(length(itoremove)>0){V<-V[-itoremove,-itoremove]}
      phylo.trimmed<-vcv2phylo(V)
      
      dat=as.matrix(dat[phylo.trimmed$tip.label,])
      LL<-mvLL(phylo.trimmed,dat,method="pic",param=list(estim=FALSE, sigma=1, check=FALSE))

    }else{
      param <- (param)
      beta<-param[seq_len(nbeta)]
      sigma<-param[nbeta+seq_len(nsigma)]
      if(!is.null(error)) errorValue <- log(param[nbeta+nsigma+1]) else errorValue <- NULL
      phylo <- BranchtransformMAPS(phylo, beta, mtot, times, fun, sigma, model, errorValue)
      phylo.trimmed<-drop.tip(phylo,phylo$tip.label[which(!phylo$tip.label%in%rownames(dat))])
      dat=as.matrix(dat[phylo.trimmed$tip.label,])

      LL<-mvLL(phylo.trimmed,dat,method="pic",param=list(estim=FALSE, sigma=1, check=FALSE))
      
    }
    if(is.na(LL$logl) | is.infinite(LL$logl)){return(1000000)}
    if(results==FALSE){
      return(-LL$logl)
    }else{
      return(list(LL=-LL$logl, mu=LL$theta, s2=sigma))
    }
  }
  
  ##------------------------------------Optimization-------------------------------##
  phyloTrans=NULL
  if(method=="BB"){
    require(BB)
    estim<-spg(par=startval,fn=function(par){clikCLIM(param=par,dat=data,phy,mtot=mtot,times=times,fun=fun,model)},control=control ,method=3, lower=lower, upper=upper)
  }else if(method=="L-BFGS-B" | method=="Nelder-Mead"){
    estim<-optim(par=startval,fn=function(par){clikCLIM(param=par,dat=data,phy,mtot=mtot,times=times,fun=fun,model)},control=control, hessian=TRUE, method=method, lower=lower, upper=upper)
  }else if(method=="fixed"){
    estim <- list()
    estim$par <- c(beta,log(sigma))
    estim$value <- clikCLIM(param=estim$par,dat=data,phy,mtot=mtot,times=times,fun=fun,model)
    estim$convergence <- 0
    phyloTrans <- BranchtransformMAPS(phy, beta, mtot, times, fun, sigma, model, errorValue=NULL)
  }
  
  # Results
  if(model=="exponential"){
    resultList<-matrix(c(estim$par[seq_len(nbeta)],exp(estim$par[nbeta+seq_len(nsigma)])),ncol=nbeta, byrow=T)
    if(constraint==FALSE) colnames(resultList)<-c(colnames(tree$mapped.edge))
    rownames(resultList)<-c("beta","sigma")
  }else{
    resultList<-matrix(estim$par[seq_len(nsigma+nbeta)],ncol=nbeta, byrow=T)
    if(constraint==FALSE) colnames(resultList)<-c(colnames(tree$mapped.edge))
    rownames(resultList)<-c("beta","sigma")
  }
  
  if(!is.null(error)) errorValue <- estim$par[nsigma+nbeta+1]^2 else errorValue <- NULL
  
  # LogLikelihood
  LL<--estim$value
  # parameter (anc + sigma)
  if(model=="exponential"){
    nparam=1+length(estim$par)
  }else{
    nparam=1+length(estim$par) #sigma estimated in optimization
  }
  
  # AIC
  AIC<--2*LL+2*nparam
  # AIC corrected
  AICc<-AIC+((2*nparam*(nparam+1))/(Ntip(phy)-nparam-1)) #Hurvich et Tsai, 1989
  
  #ancestral states estimates
  anc<-clikCLIM(param=estim$par, dat=data, phy, mtot, times, fun=fun, model, results=TRUE)$mu
  
  ##---------------------Diagnostics--------------------------------------------##
  
  if(estim$convergence==0 & diagnostic==TRUE){
    cat("\n","successful convergence of the optimizer","\n")
  }else if(estim$convergence==1 & diagnostic==TRUE){
    cat("\n","maximum limit iteration has been reached, please consider increase maxit","\n")
  }else if(diagnostic==TRUE){
    cat("\n","convergence of the optimizer has not been reached, try simpler model","\n")
  }
  
  # Hessian eigen decomposition to check the derivatives
  if(method=="BB"){
    require(numDeriv)
    #hmat<-hessian(x=estim$par, func=function(par){clikCLIM(param=par,dat=data,phy,mtot=mtot,times=times,fun=fun,model)$LL})
    #hess<-eigen(hmat)$value
    hess <- 0
  }else if(method=="L-BFGS-B" | method=="Nelder-Mead"){
    hess<-eigen(estim$hessian)$values
  }else{
    hess<-0
  }
  if(any(hess<0)){
    hess.value<-1
    if(diagnostic==TRUE){
      cat("unreliable solution has been reached, check hessian eigenvectors or try simpler model","\n")}
  }else{
    hess.value<-0
    if(diagnostic==TRUE){
      cat("a reliable solution has been reached","\n")}
  }
  
  ##-------------------Print results--------------------------------------------##
  if(echo==TRUE){
    cat("\n")
    cat("Summary results for the",model," model","\n")
    cat("LogLikelihood:","\t",LL,"\n")
    cat("AIC:","\t",AIC,"\n")
    cat("AICc:","\t",AICc,"\n")
    cat(nparam,"parameters")
    cat("\n")
    cat("Estimated rates matrix","\n")
    print(resultList)
    cat("\n")
    cat("Estimated ancestral state","\n")
    cat(anc)
    cat("\n")
    if(!is.null(error)){
      cat("\n")
      cat("Estimated error","\n")
      cat(errorValue)
      cat("\n") 
    }
  }
  
  results<-list(LogLik=LL, AIC=AIC, AICc=AICc, rates=resultList, anc=anc, convergence=estim$convergence, hess.values=hess.value, error=errorValue, param=estim$par, phyloTrans=phyloTrans)
}



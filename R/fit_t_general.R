# Fit a model for which rates depends on a time-serie curve with regime specific parameters estimates.

fit_t_general <- function(tree, data, fun, error=NULL, beta=NULL, sigma=NULL, model=c("exponential","linear"), maxdiff=0, method=c("L-BFGS-B","BB"), upper=Inf, lower=-20, control=list(maxit=20000), diagnostic=TRUE, echo=TRUE) {
  
  require(mvMORPH)
  if(!inherits(tree,"simmap")==TRUE) stop("For now only simmap-like mapped trees are allowed.","\n")
  
  # Parameters
  data<-as.matrix(data)
  method=method[1]
  rownames(data)<-tree$tip.label
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
  number_maps <- ncol(tree$mapped.edge)
  
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
  
  # Random starting value if not provided
  if(is.null(beta)){
    beta=rep(0, number_maps)
  }
  if(is.null(sigma)){
    sigma=rep(sum(pic(data,tree)^2)/n, number_maps)
  }
  
  if(model=="linear"){
    startval=c(beta,sigma)
    nbeta=length(beta)
    nsigma=length(sigma)
  }else if(model=="exponential"){
    startval=c(beta,sigma)
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
      f<-function(x, sigma, beta, funInd){sigma*exp(beta*fun[[funInd]]((mtot+maxdiff)-x))}
    }else if(model=="linear"){
      # Clim-lin function
      f<-function(x, sigma, beta, funInd){sigma+beta*fun[[funInd]]((mtot+maxdiff)-x)}
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
        bet<-beta[regimenumber]           # select the corresponding parameter for beta
        sig<-sigma[regimenumber]          # select the corresponding parameter for sigma
        bl<-currentmap[[betaval]]         # branch length under the current map
        
        tempbeta <- integrate(f, lower=age, upper=(age + bl), subdivisions=200, rel.tol = .Machine$double.eps^0.05, sigma=sig, beta=bet, funInd=regimenumber)$val
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
      sigma<-param[nbeta+seq_len(nsigma)]
      if(!is.null(error)) errorValue <- param[nbeta+nsigma+1] else errorValue <- NULL
      phylo <- BranchtransformMAPS(phylo, beta, mtot, times, fun, sigma, model, errorValue)
      
      LL<-mvLL(phylo,dat,method="pic",param=list(estim=FALSE, sigma=1, check=FALSE))
      
    }else{
      param <- (param)
      beta<-param[seq_len(nbeta)]
      sigma<-param[nbeta+seq_len(nsigma)]
      if(!is.null(error)) errorValue <- log(param[nbeta+nsigma+1]) else errorValue <- NULL
      phylo <- BranchtransformMAPS(phylo, beta, mtot, times, fun, sigma, model, errorValue)
      
      LL<-mvLL(phylo,dat,method="pic",param=list(estim=FALSE, sigma=1, check=FALSE))
      
    }
    if(is.na(LL$logl) | is.infinite(LL$logl)){return(1000000)}
    if(results==FALSE){
      return(-LL$logl)
    }else{
      return(list(LL=-LL$logl, mu=LL$theta, s2=sigma))
    }
  }
  
  ##------------------------------------Optimization-------------------------------##
  if(method=="BB"){
    require(BB)
    estim<-spg(par=startval,fn=function(par){clikCLIM(param=par,dat=data,phy,mtot=mtot,times=times,fun=fun,model)},control=control ,method=3, lower=lower, upper=upper)
  }else{
    estim<-optim(par=startval,fn=function(par){clikCLIM(param=par,dat=data,phy,mtot=mtot,times=times,fun=fun,model)},control=control, hessian=TRUE, method=method, lower=lower, upper=upper)
  }
  
  # Results
  if(model=="exponential"){
    resultList<-matrix(c(estim$par[seq_len(nsigma+nbeta)]),ncol=nbeta, byrow=T)
    colnames(resultList)<-c(colnames(tree$mapped.edge))
    rownames(resultList)<-c("beta","sigma")
  }else{
    resultList<-matrix(estim$par[seq_len(nsigma+nbeta)],ncol=nbeta, byrow=T)
    colnames(resultList)<-c(colnames(tree$mapped.edge))
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
  AICc<-AIC+((2*nparam*(nparam+1))/(n-nparam-1)) #Hurvich et Tsai, 1989
  
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
    hmat<-hessian(x=estim$par, func=clikCLIM, dat=data, phylo=phy, mtot=mtot, times=times, fun=fun, model=model)
    hess<-eigen(hmat)$value
  }else{
    hess<-eigen(estim$hessian)$values
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
  
  results<-list(LogLik=LL, AIC=AIC, AICc=AICc, rates=resultList, anc=anc, convergence=estim$convergence, hess.values=hess.value, error=errorValue, param=estim$par)
}










## ----------------- some quick and dirty tests
library(mvMORPH)

# Simulated dataset
set.seed(14)
# Generating a random tree
tree<-pbtree(n=500)

# Setting the regime states of tip species
sta<-as.vector(c(rep("Forest",250),rep("Savannah",250))); names(sta)<-tree$tip.label

# Making the simmap tree with mapped states
tree<-make.simmap(tree,sta , model="ER", nsim=1)
col<-c("blue","orange"); names(col)<-c("Forest","Savannah")

# Plot of the phylogeny for illustration
plotSimmap(tree,col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)

# Simulate the traits
sigma<-10
theta<-0
data<-mvSIM(tree, param=list(sigma=sigma, theta=theta,
                             model="BM1", nsim=1))

## Fitting the models
# BMM - Analysis with multiple rates
mvBM(tree, data)

# Make a climatic curve
require(RPANDA)
require(pspline)
data(Cetacea)
data(InfTemp)
spline_result <- sm.spline(x=InfTemp[,1],y=InfTemp[,2], df=50)
env_func <- function(t){predict(spline_result,t)}
t<-unique(InfTemp[,1])

# We build the interpolated smoothing spline function
## 1st function
env_data<-splinefun(t,env_func(t))

## 2nd function
maxtime = max(branching.times(tree))
my_fun_ebac <- function(t){
  time = (maxtime - t)
  return(time)
}

# combined
funclimeb <- list(env_data, my_fun_ebac)

# let's check
my_fun1 <- funclimeb[[1]]
curve(my_fun1, 0, maxtime)

my_fun2 <- funclimeb[[2]]
curve(my_fun2, 0, maxtime)

# Fit the models
fit_t_general(tree, data, fun=funclimeb)
fit_t_general(tree, data, fun=funclimeb, model="exponential", error=NA)
fit_t_general(tree, data, fun=funclimeb, model="exponential", error=NA, method="Nelder-Mead")
fit_t_general(tree, data, fun=funclimeb, model="linear") 
#fit_t_general(tree, data, fun=env_data, method="BB") # problem to solve



# Make a DD function?
node_Heights<-nodeHeights(tree)
maxval<- max(node_Heights)
times<-maxval-node_Heights[match(1:tree$Nnode+length(tree$tip),tree$edge[,1]),1]
# Set the root to zero
times<-maxval-times
branching<-sort(c(times,maxval))
# Extra label for the interval computation
names(branching)<-length(tree$tip.label):(length(tree$tip)*2-1)
# Number of species
tips <- length(tree$tip.label)
# Diversity index
N <- 2:tips
names(N) <- sort(times)

funN <- function(x){
  values <- as.numeric(names(N))
  res <- findInterval(x, values)
  index <- res==0
  res[index==TRUE] <- 1
  return(N[res])
}

# Let's have a look
curve(funN, 0, maxval)
  
# Make an interpolation for a smooth function? I think something along these lines...
# I think that the interpolated function will be faster than using the above one that call "findInterval"
#  t_fun <- seq(0,maxval, length.out=100)
#  funN2=splinefun(t_fun,funN(t_fun))
#  curve(funN2, 0, maxval)


# Let's try (first function for the first regime, second function for the second regime... and so on)
new_list_function <- list(Envfunction=env_data, DDfunction=funN2)

# Fit the models
fit_t_general(tree, data, fun=new_list_function)
fit_t_general(tree, data, fun=new_list_function, model="exponential", error=NA)
fit_t_general(tree, data, fun=new_list_function, model="exponential", error=NA, method="Nelder-Mead")
fit_t_general(tree, data, fun=new_list_function, model="linear") 


## NOTE: I have to change the code for the environmental dependent or the DD, because it currently assume that t=0 is the present day
# and tmax is in the past. So the funN function above is wrong... maybe a simple way to deal with that for now is to do something like that:
funN3 <- function(x) funN2(maxval-x)
curve(funN3,0,maxval)

# Not sure yet...

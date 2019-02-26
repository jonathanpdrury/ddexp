## multiple models
# mserr is expected to be the standard error of the mean for each species. Otherwise assumed to be zero.

fit_t_standard <- function(tree, data, model=c("BM1","BMM","OU1","OUM"), mserr=0, nuisance=TRUE, method="Nelder-Mead", start_val=NULL, echo=TRUE, ...){
  require(mvMORPH)
    
  model = match.arg(model)[1]
  args = list(...)
    if(is.null(args[["upper"]])) upper <- Inf else upper <- args$upper
    if(is.null(args[["lower"]])) lower <- -Inf else lower <- args$lower
    if(is.null(args[["fixedRoot"]])) fixedRoot <- TRUE else fixedRoot <- FALSE
    
  if(!inherits(tree,"simmap")==TRUE & model!="BM1" & model!="OU1") stop("For now only simmap-like mapped trees are allowed with BMM and OUM.","\n")
    else if(inherits(tree,"simmap")==TRUE) tree <- reorderSimmap(tree, order="postorder")
        else tree <- reorder(tree, order="postorder")
    
  nobs <- Ntip(tree)
  
  # Parameters
  if(!is.null(names(data))) data <- data[tree$tip.label]
    
  # Error provided?
  if(!is.null(mserr) | nuisance==TRUE){
    ## Index error
    index_error<-sapply(1:nobs, function(x){ which(tree$edge[,2]==x)})
  }
  
  # for OUM
  if(model=="OUM"){
      precalc<-mv.Precalc(tree, nb.traits=1, param=list(model="OUM", root=FALSE))
  }
  
  # for OU1 et random root
  if(model=="OU1" & fixedRoot==FALSE){
      precalc<-mv.Precalc(tree, nb.traits=1, param=list(model="OU1", root=FALSE))
  }
  
  if(model=="OUM" | model=="BMM") regimes = ncol(tree$mapped.edge)
    
  # starting values
 if(is.null(start_val)){
  sig_est <- mvLL(tree, data, param=list(estim=TRUE))$sigma # first guess
  switch(model,
        "BM1"={ 
            err <- 0.05*sig_est
            sig2 <- (1 - 0.05)*sig_est
            start_val <- log(c(sig2,err))
        },
        "BMM"={
            err <- 0.05*sig_est
            sig2 <- (1 - 0.05)*sig_est
            start_val <- log(c(rep(sig2, regimes), err)) 
        },
        "OU1"={
            times <- branching.times(tree)
            hlife <- log(2)/(max(times)/(2))
            err <- 0.05*sig_est
            sig2 <- (1 - 0.05)*sig_est
            start_val <- log(c(sig2, hlife, err))
        },
        "OUM"={
            times <- branching.times(tree)
            hlife <- log(2)/(max(times)/(2))
            err <- 0.05*sig_est
            sig2 <- (1 - 0.05)*sig_est
            start_val <- log(c(sig2, hlife, err))
        })
     if(nuisance==FALSE) start_val = start_val[-length(start_val)]
    }
    
  # log-likelihood function
    llik <- function(par, phy, data, model, error){
        switch(model,
              "BM1"={
                  sigma2 = exp(par[1])
                  phy$edge.length <- phy$edge.length*sigma2 # scaling for BM sigma
                  
                      if(!is.null(error) & nuisance==TRUE){
                          nuisance = exp(par[2])
                          phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2 + nuisance
                      }else{
                          phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2
                      }
                  
                  # ll computation
                  llik <- mvLL(phy, data, method="pic", param=list(estim=FALSE, sigma=1, check=TRUE))
              },
              "BMM"={
                 sigma2 = exp(par[1:regimes])
                 phy$edge.length <- phy$mapped.edge%*%sigma2 # scaling for BMM sigma
                  
                      if(!is.null(error) & nuisance==TRUE){
                          nuisance = exp(par[regimes+1])
                          phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2 + nuisance
                      }else{
                          phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2
                      }
                  
                  # ll computation
                  llik <- mvLL(phy, data, method="pic", param=list(estim=FALSE, sigma=1, check=TRUE)) 
              },
              "OU1"={
                  sigma2 = exp(par[1])
                  alpha = exp(par[2])
                  
                  # Tree transformation
                if(fixedRoot){ # to add later?
                    
                  phy <- phyOU(phy, alpha)
                  phy$edge.length <- phy$edge.length*sigma2 
                  
                      if(!is.null(error) & nuisance==TRUE){
                          nuisance = exp(par[3])
                          phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2 + nuisance
                      }else{
                          phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2
                      }
                  
                  # ll computation
                  llik <- mvLL(phy, data, method="pic", param=list(estim=FALSE, sigma=1, check=TRUE)) 
                  
                  }else{
                  
                   V<-.Call("mvmorph_covar_ou_random",A=precalc$C1, alpha=alpha, sigma=sigma2, PACKAGE="mvMORPH")
                    
                   if(!is.null(error) & nuisance==TRUE){
                     nuisance = exp(par[3])
                     diag(V) <- diag(V) + error^2 + nuisance
                   }else{
                     diag(V) <- diag(V) + error^2   
                   }
                  
                  # design matrix
                  W <- .Call("mvmorph_weights", nterm=as.integer(nobs), epochs=precalc$epochs, 
                           lambda=alpha, S=1, S1=1, 
                           beta=precalc$listReg, root=as.integer(0), PACKAGE="mvMORPH")
                  
                  # ll computation
                  llik <- mvLL(V, data, method="rpf", param=list(D=W))
                  
                  }
              },
              "OUM"={
                  
                  sigma2 = exp(par[1])
                  alpha = exp(par[2])
                  
                  #phy <- phyOU(phy, alpha)
                  #phy$edge.length <- phy$edge.length*sigma2 
                  #
                  #   if(!is.null(error) & nuisance==TRUE){
                  #         phy$edge.length[index_error]<-phy$edge.length[index_error]+ error^2 + nuisance
                  #    }
                  
                  # covariance matrix
                  if(fixedRoot) V<-.Call("mvmorph_covar_ou_fixed", A=precalc$C1, alpha=alpha, sigma=sigma2, PACKAGE="mvMORPH")
                   else  V<-.Call("mvmorph_covar_ou_random",A=precalc$C1, alpha=alpha, sigma=sigma2, PACKAGE="mvMORPH")
                  
                   if(!is.null(error) & nuisance==TRUE){
                     nuisance = exp(par[3])
                     diag(V) <- diag(V) + error^2 + nuisance
                   }else{
                     diag(V) <- diag(V) + error^2   
                   }
                  
                  # design matrix
                  W <- .Call("mvmorph_weights", nterm=as.integer(nobs), epochs=precalc$epochs, 
                           lambda=alpha, S=1, S1=1, 
                           beta=precalc$listReg, root=as.integer(0), PACKAGE="mvMORPH")
                  
                  # ll computation
                  llik <- mvLL(V, data, method="rpf", param=list(D=W))
              })
        
        return(llik)
    }
    
    # optimization
    if(echo==TRUE) message("Start optimization. Please wait...")
    estimModel <- optim(start_val,
                        fn = function(par) -llik(par, phy=tree, data=data, model=model, error=mserr)$logl,
                        method=method,
                        upper=upper,
                        lower=lower)
    
    # ancestral states/optimums
    if(echo==TRUE) message("Done. retrieve parameters and results...")
    theta <- as.numeric(llik(estimModel$par, phy=tree, data=data, model=model, error=mserr)$theta)
    
    # param
    param = exp(estimModel$par)
    
    if(nuisance){
    switch(model,
           "BM1"={names(param)=c("sigma2","nuisance")},
          "BMM"={names(param)=c(rep("sigma2", regimes),"nuisance")},
          "OU1"={names(param)=c("sigma2","alpha","nuisance")},
          "OUM"={names(param)=c("sigma2","alpha","nuisance")
                names(theta)=colnames(tree$mapped.edge)}) 
        }else{
    switch(model,
           "BM1"={names(param)=c("sigma2")},
          "BMM"={names(param)=rep("sigma2", regimes)},
          "OU1"={names(param)=c("sigma2","alpha")},
          "OUM"={names(param)=c("sigma2","alpha")
                names(theta)=colnames(tree$mapped.edge)})  
    }
    
    LL = -estimModel$value
    nparam = length(param)+length(theta)
    # AIC
    AIC = -2*LL+2*nparam
    # AIC corrected
    AICc = AIC+((2*nparam*(nparam+1))/(nobs-nparam-1)) 
    
    # return results
    
    results <- list(logl=LL, AIC=AIC, AICc=AICc, param=param, theta=theta, nb_param=nparam, opt=estimModel)
    
}


# function to transform the tree. We can also use rescale in geiger or other things...
 phyOU<-function(phy,alpha){
             if(alpha<=.Machine$double.eps) return(phy) # reduce to BM
             
             times <- branching.times(phy)
             names(times) <- (Ntip(phy) + 1):(Ntip(phy) + Nnode(phy))
             Tmax<-times[1]
             phy2<-phy
             
             for (i in 1:length(phy$edge.length)) {
               bl <- phy$edge.length[i]
               age <- times[which(names(times) == phy$edge[i, 1])]
               t1 <- max(times) - age
               t2 <- t1+bl
               phy2$edge.length[i] <- (1/(2*alpha))*exp(-2*alpha * (Tmax-t2)) * (1 - exp(-2 * alpha * t2)) -
                 (1/(2*alpha))*exp(-2*alpha * (Tmax-t1)) * (1 - exp(-2 * alpha * t1))
             }
             phy <- phy2
             return(phy)
           }
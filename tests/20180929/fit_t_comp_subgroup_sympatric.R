fit_t_comp_subgroup<-function(map,data,trim.class,model=c("DDexp"),par=NULL,method="Nelder-Mead",bounds=NULL,regime.map=NULL){

	if(is.null(names(data))){stop("data missing taxa names")}
	if(!is.null(dim(data))){stop("data needs to be a single trait")}
	if(is.null(bounds[["lower"]]) & is.null(bounds[["upper"]])){
        bounds$lower = -Inf
        bounds$upper = Inf
    }
    

	ClassObjectTrimmed<-CreateGeoObject_simmap.subgroup(map,trim.class)

	phylo<-ClassObjectTrimmed$map
	if(!is.null(phylo$node.label)){phylo$node.label<-NULL}
	geo.object<-ClassObjectTrimmed$geo.object	
	geo.sorted<-resortGeoObject(phylo,geo.object) 
	
	if(!is.null(regime.map)){
		#create and sort Smatrix
		regime.map.trimmed<-drop.tip.simmap(regime.map,regime.map$tip.label[which(!regime.map$tip.label%in%phylo$tip.label)])
		dist.class.object<-CreateClassObject(regime.map.trimmed)
		#this version doesn't allow for a 'both.category' (i.e., a state that is not mutually exclusive)
		dist.SMatrix<-CreateSMatrix(dist.class.object,S.cats=colnames(regime.map$mapped.edge))
		smat.sorted<-resortSMatrix(phylo,dist.SMatrix)

		#reconcile with geo object
		int<-ReconcileGeoObjectSMatrix(geo.object=geo.sorted,S.matrix=smat.sorted,phylo=phylo)
		sgeo<-int$geo.object
		smat<-int$S.matrix
	}
	
	if(length(geo.object$geography.object)<phylo$Nnode){stop("geography object cannot have more or fewer components than internode intervals in phylo")}
	if(length(data)>length(phylo$tip.label)){stop("error: some tips missing from pruned simmap")}
	if(!all(names(data)%in%phylo$tip.label)){stop("error: some tips missing from pruned simmap")}

	if(is.null(par)){par<-c(log(sqrt(var(data)/max(nodeHeights(extract.clade(phylo,getMRCA(phylo,names(data))))))),0)}
	
	if(is.null(regime.map)){
		opt<-optim(par,likelihood_subgroup_model,phylo=phylo,geography.object=geo.object,data=data,model="DDexp",method=method, lower=bounds$lower, upper=bounds$upper)
		sig2 = exp(opt$par[1])^2
		r = opt$par[2]
		z0=likelihood_subgroup_model(data=data,phylo=phylo,geography.object=geo.object,model="DDexp",par=opt$par,return.z0=TRUE)
		results<-list(model = model, LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = sig2, r = r, z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		} else {		
		##update 'par' to have enough elements for the regimes
		par=c(par[1],par[1],par[2],par[2])
		opt<-optim(par,likelihood_subgroup_model,phylo=phylo,geography.object=sgeo,data=data,r.object=smat,model="DDexp",method=method, lower=bounds$lower, upper=bounds$upper)
		sig2_1 = exp(opt$par[1])^2
		sig2_2 = exp(opt$par[2])^2
		r1 = opt$par[3]
		r2 = opt$par[4]
		z0=likelihood_subgroup_model(data=data,phylo=phylo,geography.object=geo.object,model="DDexp",par=opt$par,return.z0=TRUE)
		results<-list(model = model, LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2_1 = sig2_1, sig2_2 = sig2_2, r1 = r1, r2 = r2, z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)	
		}

}


#removing other models for testing
#	if(model=="MC"){
#		opt<-optim(par,likelihood_subgroup_model,phylo=phylo,geography.object=geo.object,data=data,model="MC",method=method, lower=bounds$lower, upper=bounds$upper)
#		sig2 = exp(opt$par[1])^2
#		S = -abs(opt$par[2])
#		z0=likelihood_subgroup_model(data=data,phylo=phylo,geography.object=geo.object,model="MC",par=opt$par,return.z0=TRUE)
#		results<-list(model = model, LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = sig2, S = S, z0 = as.numeric(z0), convergence = opt$convergence)
#		return(results)
#		}
#
#
#	if(model=="DDlin"){
#		geography.matrix<-geo.object$geography.object
#		maxN<-max(vapply(geography.matrix,function(x)max(rowSums(x)),1))
#		opt<-optim(par,likelihood_subgroup_model,phylo=phylo,geography.object=geo.object,maxN=maxN,data=data,model="DDlin",method=method, lower=bounds$lower, upper=bounds$upper)
#		sig2 = exp(opt$par[1])^2
#		b = opt$par[2]
#		z0=likelihood_subgroup_model(data=data,phylo=phylo,geography.object=geo.object,model="DDlin",par=opt$par,return.z0=TRUE,maxN=maxN)
#		results<-list(model = model, LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = sig2, b = b, z0 = as.numeric(z0), convergence = opt$convergence)
#		return(results)
#		}

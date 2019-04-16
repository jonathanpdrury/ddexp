#this is a wrapper script designed to fit  univariate EB models
# 'i.e., two-rate models vs. single-rate models (specify 'regime.map' for two-rate models, see details below)
#' models with biogeography vs. models without biogeography (specify 'geo.map' for models incorporating biogeography, see details below)
#' models with subgroup pruning vs. models without subgroup pruning (specify 'subgroup.map' for subgroup pruning algorithm, see details below)
#'
#' @param phylo: a phylogenetic tree or simmap containing all of the lineages in 'data' vector; in case of subgroup pruning, it can include other lineages
#' @param data: a named vector of continuous trait values with names corresponding to phylo$tip.labels. In case of subgroup pruning algorithm, data should only be provided for lineages that are in the subgroup at the present (i.e., length(data) won't necessarily match length(phylo$tip.label))
#' @param model: "exponential" returns exponential diversity dependent model fit, "linear" returns linear diversity dependent model fit, "both" returns a list with both fits, as elements $exponential and $linear
#' @param geo.map: an ancestral biogeography reconstruction, (matching 'phylo' object) stored as a simmap object (i.e., created using make.simmap.BGB.R)
#' @param subgroup.map: a simmap object (containing all lineages in 'phylo' object) which contains lineages in 'subgroup' specified above
#' @param subgroup: a string identifying the name of the guild ((((data just includes these? what is passed to phylo then?)))) to which to fit the model, as it appears in the simmap specified in 'subgroup.map'
#' @param regime.map: a simmap object (constructed on the exact same phylogeny passed to the 'phylo' object) which contains a reconstruction of the different rate classes (e.g., tropical and temperate lineages). Currently only tested for two-regime cases.
#' @param error: a named vector (in the same order as 'data') which specifies the standard error of trait measurements, if measurement error is to be accounted for in fit (specify as NULL if not)
#' @param beta: a vector of starting values for the slope parameter estimation
#' @param sigma: a vector of starting values for the rate parameter estimation
#' @param method: optimisation algorithm (see optim())
#' @param upper: upper bound on optimisation algorithm, if "L-BFGS-B" is chosen as 'method'
#' @param lower: lower bound on optimisation algorithm, if "L-BFGS-B" is chosen as 'method'
#' @param control: further commands passed to optim()
#' @param diagnostic: logical specifying whether diagnostic information should printed to the console
#' @param echo: logical specifying whether model fits should printed in the console

#notes for future improvement:

#--can 'phylo' be removed completely and just inferred based on maps that are passed?
#--for biogeography case, return modified results matrix

fit_t_DD<-function(phylo,data,model=c("exponential","linear","both"),regime.map=NULL,error=NULL, beta=NULL, sigma=NULL, method=c("Nelder-Mead","L-BFGS-B","BB"), upper=Inf, lower=-Inf, control=list(maxit=20000), diagnostic=FALSE, echo=FALSE){
	
	if(!model%in%c("exponential","linear","both")){ stop("model must be stated as 'exponential' , 'linear', or 'both' ")}
	
if(is.null(regime.map)){ 	# single slope version without BioGeoBEARS biogeography or subgroup pruning
	
		#check data format, names, and sorting
		
		#convert phylo to simmap object
		hold<-rep("A",length(phylo$tip.label))
		hold[1]<-"B"
		names(hold)<-phylo$tip.label
		smap<-make.simmap(phylo,hold,message=F)
		new.maps<-list()
		for(i in 1:length(phylo$edge.length)){
			new.maps[[i]]<-phylo$edge.length[i]
			names(new.maps[[i]])<-"A"
			}
		new.mapped.edge<- as.matrix(rowSums(smap$mapped.edge))
		colnames(new.mapped.edge)<-"A"	
		smap$maps<-new.maps
		smap$mapped.edge<-new.mapped.edge
	
		times = as.numeric(sort(max(branching.times(smap))-branching.times(smap)))
		
		#create function
		class.df<-return.class.df_sympatric(smap)
		new_list_function<-create.function.list(smap,times=times,df=class.df)

		#fit model
		sigma.constraint<-rep(1, dim(smap$mapped.edge)[2])
		beta.constraint<-rep(1, dim(smap$mapped.edge)[2])

		if(model%in%c("exponential","linear")){
			out<-fit_t_general(tree=smap,data=data,fun=new_list_function,class.df=class.df,input.times=times,error=error, sigma=sigma, beta=beta, model=model,method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))	
		} else {
			out.exp<-fit_t_general(tree=smap,data=data,fun=new_list_function,class.df=class.df,input.times=times,error=error, sigma=sigma, beta=beta, model="exponential",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))	
			out.lin<-fit_t_general(tree=smap,data=data,fun=new_list_function,class.df=class.df,input.times=times,error=error, sigma=sigma, beta=beta, model="linear",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))	
			out<-list(exponential.fit=out.exp,linear.fit=out.lin)
		}

		
}  else if (!is.null(regime.map)) { # two slope version without BioGeoBEARS biogeography or subgroup pruning

		#first check that regime.map and phylo and data are concordant
		if(!all(as.phylo(phylo)$tip.label == as.phylo(regime.map)$tip.label)) { stop("regime map doesn't match phylogeny")}
		if(length(data) != length(as.phylo(regime.map)$tip.label)) { stop("number of lineages in data and regime map don't match")}
		if(! all (names(data) %in% as.phylo(regime.map)$tip.label)) { stop("names of lineages in data and regime map don't match")}
		if(! all (as.phylo(regime.map)$tip.label %in% names(data)) ) { stop("names of lineages in data and regime map don't match")}
		
		class.object<-try(CreateClassObject(regime.map))
		if(class(class.object)=="try-error"){class.object<-try(CreateClassObject(regime.map,rnd=6))}
		if(class(class.object)=="try-error"){class.object<-CreateClassObject(regime.map,rnd=7)}

		class.df<-return.class.df_subgroup(regime.map,class.object)
		new_list_function<-create.function.list(regime.map,times=class.object$times,df=class.df)
				
		#fit model
		sigma.constraint<-rep(1, dim(regime.map$mapped.edge)[2])
		beta.constraint<-seq(1,by=1,length.out=dim(regime.map$mapped.edge)[2])
		
		if(model%in%c("exponential","linear")){
		out<-fit_t_general(tree=regime.map,data=data,fun=new_list_function,input.times=class.object$times,class.df=class.df,error=error, sigma=sigma, beta=beta, model=model,method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		} else{
		out.exp=fit_t_general(tree=regime.map,data=data,fun=new_list_function,input.times=class.object$times,class.df=class.df,error=error, sigma=sigma, beta=beta, model="exponential",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out.lin=fit_t_general(tree=regime.map,data=data,fun=new_list_function,input.times=class.object$times,class.df=class.df,error=error, sigma=sigma, beta=beta, model="linear",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out<-list(exponential.fit=out.exp,linear.fit=out.lin)
		}
		
}    
		return(out)
}
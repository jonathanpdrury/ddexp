#this is a wrapper script designed to fit all univariate DD models built so far
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

fit_t_DD<-function(phylo,data,model=c("exponential","linear","both"),geo.map=NULL,subgroup.map=NULL,subgroup=NULL,regime.map=NULL,error=NULL, beta=NULL, sigma=NULL, method=c("Nelder-Mead","L-BFGS-B","BB"), upper=Inf, lower=-Inf, control=list(maxit=20000), diagnostic=FALSE, echo=FALSE){
	
	if(!model%in%c("exponential","linear","both")){ stop("model must be stated as 'exponential' , 'linear', or 'both' ")}
	
if(is.null(geo.map)&&is.null(subgroup.map)&&is.null(regime.map)){ 	# single slope version without BioGeoBEARS biogeography or subgroup pruning
	
		#check data format, names, and sorting
		
#		#convert phylo to simmap object
#		hold<-rep("A",length(phylo$tip.label))
#		hold[1]<-"B"
#		names(hold)<-phylo$tip.label
#		smap<-make.simmap(phylo,hold,message=F)
#		new.maps<-list()
#		for(i in 1:length(phylo$edge.length)){
#			new.maps[[i]]<-phylo$edge.length[i]
#			names(new.maps[[i]])<-"A"
#			}
#		new.mapped.edge<- as.matrix(rowSums(smap$mapped.edge))
#		colnames(new.mapped.edge)<-"A"	
#		smap$maps<-new.maps
#		smap$mapped.edge<-new.mapped.edge
#	
#		nodeDist<-vector(mode = "numeric", length = smap$Nnode)
#  		totlen<-length(smap$tip.label)
#  		root <-totlen  + 1
#  		heights<-nodeHeights(smap)
#  		for (i in 1:dim(smap$edge)[1]){
#    		nodeDist[[smap$edge[i, 1] - totlen]] <- heights[i]
#		  }
#		nodeDist<-c(nodeDist,max(heights))
#		smap$times<-nodeDist
#		
#		#create function
#		class.df<-return.class.df_sympatric(smap)
#		new_list_function<-create.function.list_sympatric(smap,class.df)
#		
#		
#		#fit model
#		sigma.constraint<-rep(1, dim(smap$mapped.edge)[2])
#		beta.constraint<-rep(1, dim(smap$mapped.edge)[2])
#
#		out<-fit_t_general(tree=smap,data=data,fun=new_list_function,error=error, sigma=sigma, beta=beta, model=model,method=method,maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))

# THIS isn't optimising correctly for some reason; could try to reactivate if optimisation issues are workedout otherwise:

		#need to add other cases down the line
		stop("not currently implemented; please use RPANDA::fit_t_comp")

}  else if (is.null(geo.map)&&is.null(subgroup.map)&&!is.null(regime.map)) { # two slope version without BioGeoBEARS biogeography or subgroup pruning

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
		out<-fit_t_general(tree=regime.map,data=data,fun=new_list_function,class.object=class.object,class.df=class.df,error=error, sigma=sigma, beta=beta, model=model,method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		} else{
		out.exp=fit_t_general(tree=regime.map,data=data,fun=new_list_function,class.object=class.object,class.df=class.df,error=error, sigma=sigma, beta=beta, model="exponential",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out.lin=fit_t_general(tree=regime.map,data=data,fun=new_list_function,class.object=class.object,class.df=class.df,error=error, sigma=sigma, beta=beta, model="linear",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out<-list(exponential.fit=out.exp,linear.fit=out.lin)
		}
		
}  else if (is.null(subgroup.map)&&is.null(regime.map)&&!is.null(geo.map)) { # single slope version with BioGeoBEARS biogeography but no subgroup pruning

		geo.simmap<-geo.map
		hold<-CreateClassObject(geo.simmap)
		geo.class.df<-return.class.df(geo.simmap,hold)
		class.object=hold
		class.df=geo.class.df
		new_list_function <- create.function.list(geo.simmap=geo.simmap, df=class.df,times=class.object$times)

		sigma.constraint<-rep(1, dim(geo.map$mapped.edge)[2])
		beta.constraint<-rep(1, dim(geo.map$mapped.edge)[2])
		
		if(model%in%c("exponential","linear")){
		out<-fit_t_general(tree=geo.map,data=data,fun=new_list_function,class.object=class.object,class.df=class.df,error=error, sigma=sigma, beta=beta, model=model,method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		} else{
		out.exp<-fit_t_general(tree=geo.map,data=data,fun=new_list_function,class.object=class.object,class.df=class.df,error=error, sigma=sigma, beta=beta, model="exponential",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out.lin<-fit_t_general(tree=geo.map,data=data,fun=new_list_function,class.object=class.object,class.df=class.df,error=error, sigma=sigma, beta=beta, model="linear",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out<-list(exponential.fit=out.exp,linear.fit=out.lin)
		}
		
}  else if (is.null(subgroup.map)&&!is.null(regime.map)&&!is.null(geo.map)) {  # two slope version with BioGeoBEARS biogeography but no subgroup pruning
	
		#need to add other cases down the line
		stop("two-slope with biogeography not currently implemented")
	
}  else if (is.null(geo.map)&&is.null(regime.map)&&!is.null(subgroup.map)) { # single slope version with subgroup pruning but no BioGeoBEARS biogeography
		
	 	#first check  subgroup map, and phylo and data are concordant
		
		if(!all(phylo$tip.label %in% as.phylo(subgroup.map)$tip.label)) { stop("some lineages in phylogeny don't appear in subgroup map")}
		if( is.null(subgroup) || (!subgroup%in%colnames(subgroup.map$mapped.edge))){ stop("specify a subgroup that appears as a mapped regime in subgroup.map")}

		#trim subgroup simmap to species found in focal phylogeny
		subgroup.trimmed<-drop.tip.simmap(subgroup.map,subgroup.map$tip.label[which(!subgroup.map$tip.label%in%phylo$tip.label)])
		
		#trim this diet*region simmap to contain only branches in subgroup (see Drury et al. 2018 PLOS Biol for description of algorithm)
		trimclass.subgroup.trimmed<-trimSimmap(subgroup.trimmed,trim.class=subgroup)

		#create 'class object' for subsequent scripts
		class.object<-try(CreateClassObject(trimclass.subgroup.trimmed))
		
		#add some try catches to prevent crashes due to rounding issues
		if(class(class.object)=="try-error"){class.object<-try(CreateClassObject(trimclass.subgroup.trimmed,rnd=6))}
		if(class(class.object)=="try-error"){class.object<-CreateClassObject(trimclass.subgroup.trimmed,rnd=7)}

		#create class.df that has the diversity of lineages in each simmap class
		subgroup.class.df<-return.class.df_subgroup(trimclass.subgroup.trimmed,class.object)
		
		#set lineages not in subgroup to evolve via BM
		subgroup.class.df[,which(colnames(trimclass.subgroup.trimmed$mapped.edge)!=subgroup)+1]=1 
	
		#trim the subgroup.map to just lineages that are inv at the tips, adjust class.df accordingly
		
		trimclass.subgroup.trimmed.tips<-drop.tip.simmap(trimclass.subgroup.trimmed,trimclass.subgroup.trimmed$tip.label[which(!trimclass.subgroup.trimmed$tip.label%in%names(data))])
		subgroup.class.df.trimmed<-subgroup.class.df[,c(1,match(colnames(trimclass.subgroup.trimmed.tips$mapped.edge),colnames(trimclass.subgroup.trimmed$mapped.edge))+1)]		
	
		subgroup.map.region.root=max(nodeHeights(trimclass.subgroup.trimmed))
		trimclass.subgroup.trimmed.tips.root=max(nodeHeights(trimclass.subgroup.trimmed.tips))

		###	adjust the class.df and class.object if the trimmed tree has a younger root than the ancestral tree
		if(round(subgroup.map.region.root,5)!=round(trimclass.subgroup.trimmed.tips.root,5)){

			trimmed.class.object<-try(CreateClassObject(trimclass.subgroup.trimmed.tips))
			if(class(trimmed.class.object)=="try-error"){trimmed.class.object<-try(CreateClassObject(trimclass.subgroup.trimmed.tips,rnd=6))}
			if(class(trimmed.class.object)=="try-error"){trimmed.class.object<-CreateClassObject(trimclass.subgroup.trimmed.tips,rnd=7)}

			shifted.times<-trimmed.class.object$times+(subgroup.map.region.root-trimclass.subgroup.trimmed.tips.root)
			new.subgroup.class.df.trimmed<-subgroup.class.df.trimmed[c(which(round(out$times,5)==round(min(shifted.times),5)):dim(subgroup.class.df.trimmed)[1]),]
			
			new.subgroup.class.df.trimmed$interval<-c(1:dim(new.subgroup.class.df.trimmed)[1])
			class.object$times<-class.object$times[which(round(class.object$times,5)>=round(min(shifted.times),5))]-round(subgroup.map.region.root-trimclass.subgroup.trimmed.tips.root,5)
			#forces time to start at root of trimmed tree; would be better to pass times directly to new_list_function to avoid overwriting this slot of the out object, which could lead to errors
			subgroup.class.df.trimmed<-new.subgroup.class.df.trimmed
		}
	
		class.df=subgroup.class.df.trimmed
		new_list_function<-create.function.list(trimclass.subgroup.trimmed.tips,times=class.object$times,df=class.df)
				
		sigma.constraint<-rep(1, dim(trimclass.subgroup.trimmed.tips$mapped.edge)[2])
		beta.constraint<-rep(NA, dim(trimclass.subgroup.trimmed.tips$mapped.edge)[2])
		beta.constraint[which(colnames(trimclass.subgroup.trimmed.tips$mapped.edge)==subgroup)]<-1
		
		if(model%in%c("exponential","linear")){
		out<-fit_t_general(tree=trimclass.subgroup.trimmed.tips, data=data, fun=new_list_function,class.object=class.object,class.df=class.df,error=error, sigma=sigma, beta=beta, model=model,method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		} else{
		out.exp=fit_t_general(tree=trimclass.subgroup.trimmed.tips, data=data, fun=new_list_function,class.object=class.object,class.df=class.df,error=error, sigma=sigma, beta=beta, model="exponential",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out.lin=fit_t_general(tree=trimclass.subgroup.trimmed.tips, data=data, fun=new_list_function,class.object=class.object,class.df=class.df,error=error, sigma=sigma, beta=beta, model="linear",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out<-list(exponential.fit=out.exp,linear.fit=out.lin)
		}
		
}  else if (is.null(geo.map)&&!is.null(regime.map)&&!is.null(subgroup.map)) {  # two slope version with subgroup pruning but no BioGeoBEARS biogeography

		if(!all(phylo$tip.label %in% as.phylo(subgroup.map)$tip.label)) { stop("some lineages in phylogeny don't appear in subgroup map")}
		if( is.null(subgroup) || (!subgroup%in%colnames(subgroup.map$mapped.edge))){ stop("specify a subgroup that appears as a mapped regime in subgroup.map")}

		#trim subgroup.map and regime.map to species found in focal phylogeny (i.e., a phylogeny trimmed to the taxa of interest)
		
		subgroup.simmap.trimmed<-drop.tip.simmap(subgroup.map,subgroup.map$tip.label[which(!subgroup.map$tip.label%in%phylo$tip.label)])

		regime.simmap.trimmed<-drop.tip.simmap(regime.map,regime.map$tip.label[which(!regime.map$tip.label%in%phylo$tip.label)])

		##prepare simmap and fit DDexp to subgroup *regime

		#need to build a class.df where diversity represents the number of species in a regime that are ALSO in the subgroup

		class.by.class.object<-try(CreateClassbyClassObject_mvMORPH(map.guild=subgroup.simmap.trimmed,map.regime=regime.map,trim.class=subgroup))
		if(class(class.by.class.object)=="try-error"){class.by.class.object<-try(CreateClassbyClassObject_mvMORPH(map.guild=subgroup.simmap.trimmed,map.regime=regime.map,trim.class=subgroup,rnd=6))}
		if(class(class.by.class.object)=="try-error"){class.by.class.object<-CreateClassbyClassObject_mvMORPH(map.guild=subgroup.simmap.trimmed,map.regime=regime.map,trim.class=subgroup,rnd=7)}
		regime.class.df<-return.class.df_subgroup(class.by.class.object$regime.simmap,class.by.class.object$regime.class.object)

		#set lineages not in "inv" subgroup to evolve via BM
		regime.class.df[,which(colnames(class.by.class.object$regime.simmap$mapped.edge)=='Z')+1]=1

#fix naming scheme to avoid 'region', update annotation accordingly
		#trim the inv region simmap to just lineages that are inv at the tips, adjust class.df accordingly
		regime.simmap.region.trimmed<-drop.tip.simmap(class.by.class.object$regime.simmap,class.by.class.object$regime.simmap$tip.label[which(!class.by.class.object$regime.simmap$tip.label%in%names(data))])
		
		regime.class.df.trimmed<-regime.class.df[,c(1,match(colnames(regime.simmap.region.trimmed$mapped.edge),colnames(class.by.class.object$regime.simmap$mapped.edge))+1)]		

		regime.simmap.region.root=max(nodeHeights(class.by.class.object$regime.simmap))
		regime.simmap.region.trimmed.root=max(nodeHeights(regime.simmap.region.trimmed))

		###	adjust the class.df and class.object if the trimmed tree has a younger root than the ancestral tree
		if(round(regime.simmap.region.root,5)!=round(regime.simmap.region.trimmed.root,5)){
		
			trimmed.class.object<-try(CreateClassObject(regime.simmap.region.trimmed,rnd=5))
			if(class(trimmed.class.object)=="try-error"){trimmed.class.object<-try(CreateClassObject(regime.simmap.region.trimmed,rnd=6))}
			if(class(trimmed.class.object)=="try-error"){trimmed.class.object<-CreateClassObject(regime.simmap.region.trimmed,rnd=7)}

			shifted.times<-trimmed.class.object$times+(regime.simmap.region.root-regime.simmap.region.trimmed.root)
			new.regime.class.df.trimmed<-regime.class.df.trimmed[c(which(round(class.by.class.object$regime.class.object$times,5)==round(min(shifted.times),5)):dim(regime.class.df.trimmed)[1]),]
			
			new.regime.class.df.trimmed$interval<-c(1:dim(new.regime.class.df.trimmed)[1])
			class.by.class.object$regime.class.object$times<-class.by.class.object$regime.class.object$times[which(round(class.by.class.object$regime.class.object$times,5)>=round(min(shifted.times),5))]-round(regime.simmap.region.root-regime.simmap.region.trimmed.root,5)
			#forces time to start at root of trimmed tree; would be better to pass times directly to new_list_function to avoid overwriting this slot of the class.by.class.object object, which could lead to errors
			regime.class.df.trimmed<-new.regime.class.df.trimmed
		}
		class.df=regime.class.df.trimmed
		class.object=class.by.class.object$regime.class.object
		
		new_list_function<-create.function.list(regime.simmap.region.trimmed,df=class.df,times=class.object$times)
		
		sigma.constraint<-rep(1, dim(regime.simmap.region.trimmed$mapped.edge)[2])
		beta.constraint<-rep(NA, dim(regime.simmap.region.trimmed$mapped.edge)[2])
		beta.constraint[which(colnames(regime.simmap.region.trimmed$mapped.edge)!="Z")]<-1:(dim(regime.simmap.region.trimmed$mapped.edge)[2]-1)
		
		if(model%in%c("exponential","linear")){
		out<-fit_t_general(tree=regime.simmap.region.trimmed, data=data, fun=new_list_function,class.object=class.object,class.df=class.df,error=error, sigma=sigma, beta=beta, model=model,method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		} else{
		out.exp<-fit_t_general(tree=regime.simmap.region.trimmed, data=data, fun=new_list_function,class.object=class.object,class.df=class.df,error=error, sigma=sigma, beta=beta, model="exponential",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out.lin<-fit_t_general(tree=regime.simmap.region.trimmed, data=data, fun=new_list_function,class.object=class.object,class.df=class.df,error=error, sigma=sigma, beta=beta, model="linear",method=method, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		out<-list(exponential.fit=out.exp,linear.fit=out.lin)
		}
#######	
}  else if (is.null(regime.map)&&!is.null(subgroup.map)&&!is.null(geo.map)) { # single slope version with subgroup pruning and BioGeoBEARS biogeography
	
		#need to add other cases down the line (check 20180802_testing.R and CreateGeoByClassObject.R for inspiration)
		stop("single slope version with subgroup pruning and BioGeoBEARS biogeography not yet implemented")
	
	}  else if (!is.null(regime.map)&&!is.null(subgroup.map)&&!is.null(geo.map)) {  # two slope version with subgroup pruning and BioGeoBEARS biogeography
	
		#need to add other cases down the line
		stop("two slope version with subgroup pruning and BioGeoBEARS biogeography not yet implemented")
	}
		return(out)
}
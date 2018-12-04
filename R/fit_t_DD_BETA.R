#this is a wrapper script designed to fit all univariate DD models built so far
#i.e., two rate models vs. single rate models
#models with biogeography vs. models without biogeography
#models with subgroup pruning vs. models without subgroup pruning


source('fit_t_general.R')
source('generalized_functions.R')

fit_t_DD<-function(phylo,data,model=c("exponential","linear"),geo.map=NULL,subgroup=NULL,subgroup.map=NULL,error=NULL, beta=NULL, sigma=NULL, method=c("L-BFGS-B","BB"), upper=Inf, lower=-20, control=list(maxit=20000), diagnostic=FALSE, echo=FALSE){
	
	if(!model%in%c("exponential","linear")){ stop("model must be stated as 'exponential' or 'linear'")}
	
	if(is.null(geo.map)&&is.null(subgroup.map)){
	
		#check data format, names, and sorting
		
		#convert phylo to simmap object
		hold<-rep("A",length(phylo)$tip.label)
		hold[1]<-"B"
		names(hold)<-phylo$tip.label
		smap<-make.simmap(phylo,hold)
		new.maps<-list()
		for(i in 1:length(phylo$edge.length)){
			new.maps[[i]]<-phylo$edge.length[i]
			names(new.maps[[i]])<-"A"
			}
		new.mapped.edge<- as.matrix(rowSums(smap$mapped.edge))
		colnames(new.mapped.edge)<-"A"	
		smap$maps<-new.maps
		smap$mapped.edge<-new.mapped.edge
	
		nodeDist<-vector(mode = "numeric", length = smap$Nnode)
  		totlen<-length(smap$tip.label)
  		root <-totlen  + 1
  		heights<-nodeHeights(smap)
  		for (i in 1:dim(smap$edge)[1]){
    		nodeDist[[smap$edge[i, 1] - totlen]] <- heights[i]
		  }
		nodeDist<-c(nodeDist,max(heights))
		smap$times<-nodeDist
		
		#create function
		class.df<-return.class.df_sympatric(smap)
		new_list_function<-create.function.list_sympatric(smap,class.df)
		
		#calculate maxN if DDlin, set to NA if DDexp
		
		maxN<-ifelse(model=="linear",max(class.df[,2]),NA)
		
		#fit model
		
		out<-fit_t_general(tree=smap,data=data,fun=new_list_function,error=error, sigma=sigma, beta=beta, model=model,method=method,maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=c(1,1), beta=c(1,1)))

	} else{ #need to add other cases down the line
	
		stop("only sympatric clades without subgroup pruning possible in current wrapper script")
	
	}  

		return(out)
}
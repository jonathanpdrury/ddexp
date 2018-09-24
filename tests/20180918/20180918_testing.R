require(deSolve)
require(geiger)
require(phytools)
library(corpcor)
require(methods)
require(tools)
require(RPANDA)
require(data.table)
require(mvMORPH)

source('fit_t_general.R')
source('generalized_functions.R')
source('CreateGeobyClassObject_mvMORPH.R')
source('make.simmap.BGB.R')

source('~/Dropbox/Scripts/R scripts/stratified_BGB_to_tables.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/trimSimmap.R', chdir = TRUE)

load('~/ddexp/data/motmot.data.Rd')

smap<-stratified_BGB_to_tables(motmot.data[[1]],motmot.data[[2]],1)
geo.object<-CreateBioGeoB_Object_subclade(anc.phylo=motmot.data[[1]],subclade.phylo=motmot.data[[1]],ana.events=smap$ana.int,clado.events=smap$clado.int)


##create class simmap where root will conflict

trait<-c(rep("A",5),rep("B",5))
names(trait)<-motmot.data[[1]]$tip.label
#smap1<-make.simmap(motmot.data[[1]],trait,model="ARD")
#
#smap1$maps[[15]]<-c(1,12, 13.93837-13)
#names(smap1$maps[[15]])<-c("B","A","B")
#
#save(smap1,file="smap1.RData")
load('~/ddexp/tests/20180918/smap1.RData')

smap<-stratified_BGB_to_tables(motmot.data[[1]],motmot.data[[2]],1)
out<-CreateGeobyClassObject_mvMORPH(phylo=motmot.data[[1]],map=smap1,ana.events=smap$ana.int,clado.events=smap$clado.int,trim.class="A",rnd=5)
	
geo.class.df<-return.class.df(out$geo.simmap,out$geo.class.object)
geo.class.df[,which(colnames(out$geo.simmap$mapped.edge)=='Z')+1]=1 #this is the change from the 0802 tests

M<-motmot.data[[3]][which(trait=="A")]
geo.simmap.trimmed<-drop.tip.simmap(out$geo.simmap,out$geo.simmap$tip.label[which(!out$geo.simmap$tip.label%in%names(M))])
geo.class.df.trimmed<-geo.class.df[,c(1,match(colnames(geo.simmap.trimmed$mapped.edge),colnames(out$geo.simmap$mapped.edge))+1)]		

geo.simmap.root=max(nodeHeights(out$geo.simmap))
geo.simmap.trimmed.root=max(nodeHeights(geo.simmap.trimmed))

##roots aren't the same age (as desired); now designing a workaround

if(round(geo.simmap.root,5)!=round(geo.simmap.trimmed.root,5)){

	trimmed.geoclass.object<-CreateClassObject(geo.simmap.trimmed)
	shifted.times<-trimmed.geoclass.object$times+(geo.simmap.root-geo.simmap.trimmed.root)
	new.geo.class.df.trimmed<-geo.class.df.trimmed[c(which(round(out$geo.class.object$times,5)==round(min(shifted.times),5)):dim(geo.class.df.trimmed)[1]),]
	new.geo.class.df.trimmed$interval<-c(1:dim(new.geo.class.df.trimmed)[1])
	out$geo.class.object$times<-out$geo.class.object$times[which(round(out$geo.class.object$times,5)>=round(min(shifted.times),5))]-round(geo.simmap.root-geo.simmap.trimmed.root,5)
	#forces time to start at root of trimmed tree; would be better to pass times directly to new_list_function to avoid overwriting this slot of the out object, which could lead to errors
}
	
new_list_function<-create.function.list(geo.simmap.trimmed,out$geo.class.object,new.geo.class.df.trimmed)
o5<-fit_t_general(geo.simmap.trimmed, M, fun=new_list_function,diagnostic=T,echo=T)


# successful convergence of the optimizer 
#a reliable solution has been reached 
#
#Summary results for the exponential  model 
#LogLikelihood: 	 -2.54914 
#AIC: 	 11.09828 
#AICc: 	 35.09828 
#3 parameters
#Estimated rates matrix 
#            [,1]
#beta  0.03433709
#sigma 0.01042773
#
#Estimated ancestral state 
#4.502268



##now check against RPANDA (which does more straightforward manipulation of vcv matrix)

fit_t_comp_subgroup(full.phylo=motmot.data[[1]],ana.events=smap$ana.int,clado.events=smap$clado.int,map=smap1,data=M,trim.class="A",model="DDexp",par=NULL,method="Nelder-Mead",bounds=NULL)

#fit_t_comp_subgroup(full.phylo=motmot.data[[1]],ana.events=smap$ana.int,clado.events=smap$clado.int,map=smap1,data=M,trim.class="A",model="DDexp",par=NULL,method="Nelder-Mead",bounds=NULL)
#$model
#[1] "DDexp"
#
#$LH
#[1] -2.978867
#
#$aic
#[1] 11.95773
#
#$aicc
#[1] 23.95773
#
#$free.parameters
#[1] 3
#
#$sig2
#[1] 0.005142511
#
#$r
#[1] 0.2858203
#
#$z0
#[1] 4.492169
#
#$convergence
#[1] 0


## not working yet


##I wonder if this could be a result of the effect of trimming off the root--in other words the variance-covariance structure of the data changes if the branch extending from the root is still present versus trimmed

##Not sure how to test this with the current approach--could perhaps try to analytically calculate vcv under different approaches and see what is returned by function?
##With old approach, can check whether VCV matrix returned with slope = 0 is aligned with vcv for trimmed tree or for subsection of hte full tree--can also look to see whether LH returned from fitContinuous() is the same under likelihood_t_subgroup

full.phylo=motmot.data[[1]]
ana.events=smap$ana.int
clado.events=smap$clado.int
map=smap1
data=M
trim.class="A"
model="DDexp"
par=NULL
method="Nelder-Mead"
bounds=NULL
stratified=FALSE

source('~/Dropbox/Scripts/R scripts/resortGeoObject.R', chdir = TRUE)

	if(is.null(names(data))){stop("data missing taxa names")}
	if(!is.null(dim(data))){stop("data needs to be a single trait")}
	is_tip <- full.phylo$edge[,2] <= length(full.phylo$tip.label)
	if(sum(diff(full.phylo$edge[is_tip, 2])<0)>0){ stop('fit_t_comp_subgroup cannot be used with ladderized full.phylogenies')}
	
	if(is.null(bounds[["lower"]]) & is.null(bounds[["upper"]])){
        bounds$lower = -Inf
        bounds$upper = Inf
    }
    
	GeoByClassObject<-CreateGeobyClassObject(full.phylo,map,trim.class,ana.events,clado.events,stratified=stratified)

	phylo<-GeoByClassObject$map
	if(!is.null(phylo$node.label)){phylo$node.label<-NULL}
	geo.object<-GeoByClassObject$geo.object
	#geo.sorted<-.resortGeoObject(phylo,geo.object) 

	if(length(geo.object$geography.object)<phylo$Nnode){stop("geography object cannot have more or fewer components than internode intervals in phylo")}
	if(length(data)>length(phylo$tip.label)){stop("error: some tips missing from pruned simmap")}
	if(!all(names(data)%in%phylo$tip.label)){stop("error: some tips missing from pruned simmap")}

geo.sorted<-resortGeoObject(phylo,geography.object) 

params<-c(0,log(1),0)
ddexp.ob<-createModel_DDexp_geo(phylo,geo.sorted)
tipdistribution <- getTipDistribution(ddexp.ob, c(params))            
V_RPANDA<-tipdistribution$Sigma
V_PHYLO<-vcv.phylo(phylo)

#so the VCV matrices under the *full* phylo match what occurs under RPANDAs trimming algorithm
#which is to say, does NOT match vcv matrix under simple BM model, which would instead be

V_PHYLO.trimmed<-vcv.phylo(drop.tip(phylo,"Electron_platyrhynchum"))

##So the BM ML estimate should return a different likelihood under subgroup model

bm.mle<-fitContinuous(drop.tip(phylo,"Electron_platyrhynchum"),M)

#GEIGER-fitted comparative model of continuous data
# fitted ‘BM’ model parameters:
#	sigsq = 0.011486
#	z0 = 4.504002
#
# model summary:
#	log-likelihood = -2.549870
#	AIC = 9.099740
#	AICc = 15.099740
#	free parameters = 2
#
#Convergence diagnostics:
#	optimization iterations = 100
#	failed iterations = 0
#	frequency of best fit = 1.00
#
# object summary:
#	'lik' -- likelihood function
#	'bnd' -- bounds for likelihood search
#	'res' -- optimization iteration summary
#	'opt' -- maximum likelihood parameter estimates

par<-c(log(sqrt(0.011486)),0)
likelihood_subgroup_model(data=M,phylo=phylo,geography.object=geo.object,model="DDexp",par=par)

#         [,1]
#[1,] 3.023324
#attr(,"logarithm")
#[1] TRUE


#whereas the approach above (the mvMORPH approach) should return the right likelihood
##(try this)

# Use the function in fit_t_general to recover the transformed tree (i.e. we don't optimize just fix the parameters)

test<-fit_t_general(geo.simmap.trimmed, M, fun=new_list_function,diagnostic=T,echo=T,method="fixed")
tempTree <- test$phyloTrans # retrieve the transformed tree
class(tempTree) <- "phylo"  # I change it to a phylo object because the simmap structure is affected...
plot(tempTree) # Let's plot it

#
#Summary results for the exponential  model 
#LogLikelihood: 	 -2.54987 
#AIC: 	 11.09974 
#AICc: 	 35.09974 
#3 parameters
#Estimated rates matrix 
#            [,1]
#beta  0.00000000
#sigma 0.01148575


#> likelihood_t_DD(drop.tip(phylo,"Electron_platyrhynchum"),M,par=c(log(0.01148575),0),model="DDexp")
#[1] 2.54987


##SO, the mvMORPH approach is returning the expected likelihood when sig2 == 0; I think that it is likely returning reliable estimates in this scenario. Are there other tests to run? Calculate expected VCV matrix by hand for this simple scenario.


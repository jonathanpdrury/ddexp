##this script checks wither subgroup pruning algorithm using mvMORPH approach is working

source('~/Dropbox/Scripts/R scripts/trimSimmap.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/CreateClassObject.R', chdir = TRUE)
source('~/ddexp/tests/20180929/fit_t_comp_subgroup_sympatric.R', chdir = TRUE)
source('~/ddexp/tests/20180929/likelihood_subgroup_model_sympatric.R', chdir = TRUE)
source('~/ddexp/R/CreateBioGeoB_Object_subclade.R', chdir = TRUE)
source('~/ddexp/R/CreateClassbyClassObject_mvMORPH.r', chdir = TRUE)
source('~/ddexp/R/CreateGeobyClassObject_mvMORPH.R', chdir = TRUE)
source('~/ddexp/R/fit_t_general_subgroup.R', chdir = TRUE)
source('~/ddexp/R/fit_t_general.R', chdir = TRUE)
source('~/ddexp/R/generalized_functions.R', chdir = TRUE)
source('~/ddexp/R/make.simmap.BGB.R', chdir = TRUE)
source('~/ddexp/R/stratified_BGB_to_tables.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/DDexp_geo_ADiag.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/resortGeoObject.R', chdir = TRUE)
source('~/ddexp/tests/20180929/CreateGeoObject_simmap.subgroup.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/CreateSMatrix.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/ReconcileGeoObjectSMatrix.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/resortSMatrix.R', chdir = TRUE)
source('~/Dropbox/VCV rescale improving (MC)/DDexp multi-slope model/DDexpMulti_geo_ADiag.R', chdir = TRUE)
source('~/ddexp/tests/20180929/createDDM_BETA.R', chdir = TRUE)

load('~/Dropbox/Manuscripts/PLOS Biology Tanager Project/data/stochastic maps for MCC tree/diet.simmaps.RData')

master<-read.csv("~/Dropbox/Manuscripts/PLOS Biology Tanager Project/data/MASTER_tandata.csv")
smap<-diet.simmaps[[1]]
smap$node.label<-NULL
	
fruit<-subset(master,diet=="fruit")

M<-fruit$ln.mass
names(M)<-fruit$treename
M<-subset(M,M!="NA" & names(M) %in% smap$tip.label)

#make tropical/temperate dummy simmap

tree<-as.phylo(smap)

trait<-rep(c(rep("tropical",20),rep("temperate",12)),40)[1:351]
names(trait)<-tree$tip.label
set.seed(123)
smap2<-make.simmap(tree,trait,model="ARD")
save(smap2,file="smap2.RData")

##fitting with RPANDA-esque approach:

o1<-fit_t_comp_subgroup(smap,M,trim.class="fruit",regime.map=smap2)

#> o1
#$model
#[1] "DDexp"
#
#$LH
#[1] -16.68427
#
#$aic
#[1] 39.36854
#
#$aicc
#[1] 39.60155
#
#$free.parameters
#[1] 3
#
#$sig2_1
#[1] 1.459219e-14
#
#$sig2_2
#[1] 0.02209943
#
#$r1
#[1] 1.188298
#
#$r2
#[1] -0.03965939
#
#$z0
#[1] 3.215475
#
#$convergence
#[1] 0


##fitting with mvMORPH-esque approach:

#need to build a class.df where diversity represents the number of species in a REGION that are ALSO in the fruit category
		
out<-CreateClassbyClassObject_mvMORPH(map.guild=smap,map.regime=smap2,trim.class="fruit")
regime.class.df<-return.class.df(out$regime.simmap,out$regime.class.object)

#set lineages not in "fruit" subgroup to evolve via BM
regime.class.df[,which(colnames(out$regime.simmap$mapped.edge)=='Z')+1]=1

#trim the fruit region simmap to just lineages that are fruit at the tips, adjust class.df accordingly
regime.simmap.region.trimmed<-drop.tip.simmap(out$regime.simmap,out$regime.simmap$tip.label[which(!out$regime.simmap$tip.label%in%names(M))])
regime.class.df.trimmed<-regime.class.df[,c(1,match(colnames(regime.simmap.region.trimmed$mapped.edge),colnames(out$regime.simmap$mapped.edge))+1)]		

regime.simmap.region.root=max(nodeHeights(out$regime.simmap))
regime.simmap.region.trimmed.root=max(nodeHeights(regime.simmap.region.trimmed))

####	adjust the class.df and class.object if the trimmed tree has a younger root than the ancestral tree
#if(round(regime.simmap.region.root,5)!=round(regime.simmap.region.trimmed.root,5)){
#
#	trimmed.class.object<-CreateClassObject(regime.simmap.region.trimmed)
#	shifted.times<-trimmed.class.object$times+(regime.simmap.region.root-regime.simmap.region.trimmed.root)
#	new.regime.class.df.trimmed<-regime.class.df.trimmed[c(which(round(out$regime.class.object$times,5)==round(min(shifted.times),5)):dim(regime.class.df.trimmed)[1]),]
#	
#	new.regime.class.df.trimmed$interval<-c(1:dim(new.regime.class.df.trimmed)[1])
#	out$regime.class.object$times<-out$regime.class.object$times[which(round(out$regime.class.object$times,5)>=round(min(shifted.times),5))]-round(regime.simmap.region.root-regime.simmap.region.trimmed.root,5)
#	#forces time to start at root of trimmed tree; would be better to pass times directly to new_list_function to avoid overwriting this slot of the out object, which could lead to errors
#}

new_list_function<-create.function.list(regime.simmap.region.trimmed,out$regime.class.object,regime.class.df.trimmed)
o3B<-fit_t_general(regime.simmap.region.trimmed, M, fun=new_list_function,diagnostic=T,echo=T,constraint=FALSE)			 

#successful convergence of the optimizer 
#a reliable solution has been reached 
#
#Summary results for the exponential  model 
#LogLikelihood: 	 -18.14244 
#AIC: 	 50.28489 
#AICc: 	 52.43873 
#7 parameters
#Estimated rates matrix 
#         temperate          Z    tropical
#beta  0.0452046634 0.68738265 -0.02179255
#sigma 0.0004205114 0.01206599  0.01463595

#One issue is that this script estimates parameters for the "Z" part of the map (where lineages are in allopatry)
#Need to figure out proper workaround

fit_t_general(regime.simmap.region.trimmed, beta=c(1.188298,0,-0.03965939),sigma=c(1.459219e-14,0.02209943,0.02209943), M, fun=new_list_function,diagnostic=T,echo=T,constraint=FALSE,method="fixed")


###getting *very* different likelihoods--need to figure out where the issue is
###also need to think about how to deal with "Z" parameter estimation
##>>>could edit fit_t_general to force beta for Z to be 0; this is accomplished by setting regime.class.df column to 0
##>>>maybe also to force sig2 to be equal across regimes?

##ultimately, this approach is different from RPANDA approach; would need to conduct a wholescale rewrite to conduct proper comparison tests;
##currently, two regime approach works fine, as does subgroup pruning, so I think it is safe to go forward with this approach; 
##perhaps running a small simulation study to test the ability for hte method to detect things?



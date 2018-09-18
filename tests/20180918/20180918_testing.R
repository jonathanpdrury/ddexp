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

load('~/ddexp/data/motmot.data.Rd')

smap<-stratified_BGB_to_tables(motmot.data[[1]],motmot.data[[2]],1)
geo.object<-CreateBioGeoB_Object_subclade(anc.phylo=motmot.data[[1]],subclade.phylo=motmot.data[[1]],ana.events=smap$ana.int,clado.events=smap$clado.int)


##create class simmap where root will conflict

trait<-c(rep("A",5),rep("B",5))
names(trait)<-motmot.data[[1]]$tip.label
smap1<-make.simmap(motmot.data[[1]],trait,model="ARD")

smap1$maps[[15]]<-c(1,12, 13.93837-13)
names(smap1$maps[[15]])<-c("B","A","B")

save(smap1,file="smap1.RData")

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
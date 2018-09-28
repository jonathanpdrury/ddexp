##this script checks wither subgroup pruning algorithm using mvMORPH approach is working

source('~/Dropbox/Scripts/R scripts/trimSimmap.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/CreateClassObject.R', chdir = TRUE)
source('~/ddexp/tests/20180928/fit_t_comp_subgroup_sympatric.R', chdir = TRUE)
source('~/ddexp/tests/20180928/likelihood_subgroup_model_sympatric.R', chdir = TRUE)
source('~/ddexp/R/CreateBioGeoB_Object_subclade.R', chdir = TRUE)
source('~/ddexp/R/CreateClassbyClassObject_mvMORPH.r', chdir = TRUE)
source('~/ddexp/R/CreateGeobyClassObject_mvMORPH.R', chdir = TRUE)
source('~/ddexp/R/fit_t_general_subgroup.R', chdir = TRUE)
source('~/ddexp/R/fit_t_general.R', chdir = TRUE)
source('~/ddexp/R/generalized_functions.R', chdir = TRUE)
source('~/ddexp/R/make.simmap.BGB.R', chdir = TRUE)
source('~/ddexp/R/stratified_BGB_to_tables.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/DDlin_geo_ADiag.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/MC_geo_PM.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/resortGeoObject.R', chdir = TRUE)
source('~/ddexp/tests/20180928/CreateGeoObject_simmap.subgroup.R', chdir = TRUE)

load('~/Dropbox/Manuscripts/PLOS Biology Tanager Project/data/stochastic maps for MCC tree/diet.simmaps.RData')

master<-read.csv("~/Dropbox/Manuscripts/PLOS Biology Tanager Project/data/MASTER_tandata.csv")
smap<-diet.simmaps[[1]]
smap$node.label<-NULL
	
fruit<-subset(master,diet=="fruit")

M<-fruit$ln.mass
names(M)<-fruit$treename
M<-subset(M,M!="NA" & names(M) %in% smap$tip.label)

#fit with RPANDA-esque approach

#> o1
#$model
#[1] "DDexp"
#
#$LH
#[1] -19.63744
#
#$aic
#[1] 45.27489
#
#$aicc
#[1] 45.5079
#
#$free.parameters
#[1] 3
#
#$sig2
#[1] 0.01107849
#
#$r
#[1] -0.01413724
#
#$z0
#[1] 3.102948
#
#$convergence
#[1] 0



## NOW using mvMORPH-esque approach

		##prepare simmap and fit DDexp to subgroup

		fruit.simmap<-trimSimmap(smap,trim.class="fruit")
		
		#create 'class object' for subsequent scripts
		out<-CreateClassObject(fruit.simmap)
		
		#create class.df that has the diversity of lineages in each simmap class
		fruit.class.df<-return.class.df_subgroup(fruit.simmap,out)
		
		#set lineages not in "fruit" subgroup to evolve via BM
		fruit.class.df[,which(colnames(fruit.simmap$mapped.edge)!='fruit')+1]=1 
	
		#trim the fruit region simmap to just lineages that are fruit at the tips, adjust class.df accordingly
		fruit.simmap.trimmed<-drop.tip.simmap(fruit.simmap,fruit.simmap$tip.label[which(!fruit.simmap$tip.label%in%names(M))])
		fruit.class.df.trimmed<-fruit.class.df[,c(1,match(colnames(fruit.simmap.trimmed$mapped.edge),colnames(fruit.simmap$mapped.edge))+1)]		

		fruit.simmap.root=max(nodeHeights(fruit.simmap))
		fruit.simmap.trimmed.root=max(nodeHeights(fruit.simmap.trimmed))

#		###	adjust the class.df and class.object if the trimmed tree has a younger root than the ancestral tree
#		if(round(fruit.simmap.root,5)!=round(fruit.simmap.trimmed.root,5)){
#
#			trimmed.class.object<-CreateClassObject(fruit.simmap.region.trimmed)
#			shifted.times<-trimmed.class.object$times+(fruit.simmap.region.root-fruit.simmap.region.trimmed.root)
#			new.fruit.class.df.trimmed<-fruit.class.df.trimmed[c(which(round(out$times,5)==round(min(shifted.times),5)):dim(fruit.class.df.trimmed)[1]),]
#			
#			new.fruit.class.df.trimmed$interval<-c(1:dim(new.fruit.class.df.trimmed)[1])
#			out$times<-out$times[which(round(out$times,5)>=round(min(shifted.times),5))]-round(fruit.simmap.region.root-fruit.simmap.region.trimmed.root,5)
#			#forces time to start at root of trimmed tree; would be better to pass times directly to new_list_function to avoid overwriting this slot of the out object, which could lead to errors
#		}
		
		new_list_function<-create.function.list(fruit.simmap.trimmed,out,fruit.class.df.trimmed)
		
		o3<-fit_t_general(fruit.simmap.trimmed, M, fun=new_list_function,diagnostic=T,echo=T)			 

#successful convergence of the optimizer 
#a reliable solution has been reached 
#
#Summary results for the exponential  model 
#LogLikelihood: 	 -19.64202 
#AIC: 	 45.28404 
#AICc: 	 45.71261 
#3 parameters
#Estimated rates matrix 
#             [,1]
#beta  -0.01406261
#sigma  0.01104329
#
#Estimated ancestral state 
#3.103023


##Looks good!
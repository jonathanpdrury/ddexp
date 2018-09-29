##This script fits BM, OU, DD, plus same models with latitude as categorical variable
##while also doing subgroup trimming

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
source('trimSimmap.R')
source('CreateClassObject.R')


master<-read.csv("BIRDS_MASTER.csv")
tree<-read.tree("treefile_allbirds_simple.tree")

load('diet.simmap.allbirds.RData')
load('all.simmap.trop.RData')


#diet.simmaps were generated with all data, whereas biogeography was limited to a subset of species; need to trim to reflect this
diet.simmap.allbirds.i<-diet.simmap.allbirds[[1]]
diet.simmap.allbirds.i<-drop.tip.simmap(diet.simmap.allbirds.i,diet.simmap.allbirds.i$tip.label[which(!diet.simmap.allbirds.i$tip.label%in%tree$tip.label)])

trop.simmap.allbirds.i<-all.simmap.trop[[1]]
trop.simmap.allbirds.i<-drop.tip.simmap(trop.simmap.allbirds.i,trop.simmap.allbirds.i$tip.label[which(!trop.simmap.allbirds.i$tip.label%in%tree$tip.label)])

region.vec<-c("AB","CD","E","F","G","H","I","J","K","L")

res.mat<-matrix(nrow=7,ncol=21)
colnames(res.mat)<-c("region","trait","diet.category","no","BM1.lnL","BM1.AIC","BM1.sig","BM1.z0","BM1.convergence","OU1.lnL","OU1.AIC","OU1.sig","OU1.alpha","OU1.z0","OU1.convergence","geo.map","DDexp.lnL","DDexp.sig","DDexp.r","DDexp.z0","DDexp.convergence")
outlist=list()

fruit<-subset(master,diet=="fruit")


for(i in 1:length(region.vec)){

	#trim tree down to region of interest
	
	if(i==1){
	region<-subset(master,A==1 | B==1)
	region.tree<-drop.tip(tree,tree$tip.label[which(!tree$tip.label%in%region$Species_name_from_birdtree)])
	}
	if(i==2){
	region<-subset(master,C==1 | D==1)
	region.tree<-drop.tip(tree,tree$tip.label[which(!tree$tip.label%in%region$Species_name_from_birdtree)])
	}
	if(i>2){
	eval(parse(text=paste0("region<-subset(master,",region.vec[i],"==1)")))
	region.tree<-drop.tip(tree,tree$tip.label[which(!tree$tip.label%in%region$Species_name_from_birdtree)])
	}
	
	##need to further trim until simmap is constructed for expanded set of species
	fruit<-fruit[fruit$Species_name_from_birdtree%in%diet.simmap.allbirds.i$tip.label,]
	
	for(j in 1:7){
		M<-fruit[,90+j]
		names(M)<-fruit[,53]
		M<-subset(M,M!="NA" & names(M) %in% region.tree$tip.label)
		nc<-name.check(region.tree,M)
		nc
		if(is.list(nc)){
			subtree<-drop.tip(region.tree,nc$tree_not_data)
			} else{
			subtree<-region.tree
			}
		M<-M[subtree$tip.label]
		
		##fit BM and OU first
		o1<-mvBM(subtree,M,model=c("BM1"),optimization="subplex")
		o2<-mvOU(subtree,M,model=c("OU1"),method="sparse", optimization="subplex")
		
		##prepare simmap and fit DDexp to subgroup
		
		#trim diet simmap to species found in focal region
		diet.simmap.region<-drop.tip.simmap(diet.simmap.allbirds.i,diet.simmap.allbirds.i$tip.label[which(!diet.simmap.allbirds.i$tip.label%in%region.tree$tip.label)])
		
		#trim this diet*region simmap to contain only branches in "fruit" subgroup (see Drury et al. 2018 PLOS Biol for description of algorithm)
		fruit.simmap.region<-trimSimmap(diet.simmap.region,trim.class="fruit")
		
		#create 'class object' for subsequent scripts
		out<-CreateClassObject(fruit.simmap.region)
		
		#create class.df that has the diversity of lineages in each simmap class
		fruit.class.df<-return.class.df_subgroup(fruit.simmap.region,out)
		
		#set lineages not in "fruit" subgroup to evolve via BM
		fruit.class.df[,which(colnames(fruit.simmap.region$mapped.edge)!='fruit')+1]=1 
	
		#trim the fruit region simmap to just lineages that are fruit at the tips, adjust class.df accordingly
		fruit.simmap.region.trimmed<-drop.tip.simmap(fruit.simmap.region,fruit.simmap.region$tip.label[which(!fruit.simmap.region$tip.label%in%names(M))])
		fruit.class.df.trimmed<-fruit.class.df[,c(1,match(colnames(fruit.simmap.region.trimmed$mapped.edge),colnames(fruit.simmap.region$mapped.edge))+1)]		
	
		fruit.simmap.region.root=max(nodeHeights(fruit.simmap.region))
		fruit.simmap.region.trimmed.root=max(nodeHeights(fruit.simmap.region.trimmed))

		###	adjust the class.df and class.object if the trimmed tree has a younger root than the ancestral tree
		if(round(fruit.simmap.region.root,5)!=round(fruit.simmap.region.trimmed.root,5)){

			trimmed.class.object<-CreateClassObject(fruit.simmap.region.trimmed)
			shifted.times<-trimmed.class.object$times+(fruit.simmap.region.root-fruit.simmap.region.trimmed.root)
			new.fruit.class.df.trimmed<-fruit.class.df.trimmed[c(which(round(out$times,5)==round(min(shifted.times),5)):dim(fruit.class.df.trimmed)[1]),]
			
			new.fruit.class.df.trimmed$interval<-c(1:dim(new.fruit.class.df.trimmed)[1])
			out$times<-out$times[which(round(out$times,5)>=round(min(shifted.times),5))]-round(fruit.simmap.region.root-fruit.simmap.region.trimmed.root,5)
			#forces time to start at root of trimmed tree; would be better to pass times directly to new_list_function to avoid overwriting this slot of the out object, which could lead to errors
		}
		
		new_list_function<-create.function.list(fruit.simmap.region.trimmed,out,fruit.class.df.trimmed)
		o3<-fit_t_general(fruit.simmap.region.trimmed, M, fun=new_list_function,diagnostic=T,echo=T)			 
		DDexpgeo_mvMORPH.lnL<-o3$LogLik
		DDexpgeo_mvMORPH.sig2<-o3$rates[2,1]
		DDexpgeo_mvMORPH.r<-o3$rates[1,1]
		DDexpgeo_mvMORPH.z0<-o3$anc
		DDexpgeo_mvMORPH.conv<-o3$convergence
		
		### NOW fitting models that vary according to tropical/temperate membership
		
		trop.simmap.trimmed<-drop.tip.simmap(trop.simmap.allbirds.i,trop.simmap.allbirds.i$tip.label[which(!trop.simmap.allbirds.i$tip.label%in%subtree$tip.label)])
		o1B<-mvBM(trop.simmap.trimmed,M,model=c("BMM"),optimization="subplex")
		o2B<-mvOU(trop.simmap.trimmed,M,model=c("OUM"),optimization="subplex")

		##prepare simmap and fit DDexp to subgroup *regime

		#need to build a class.df where diversity represents the number of species in a REGION that are ALSO in the fruit category
		
		out<-CreateClassbyClassObject_mvMORPH(map.guild=diet.simmap.region,map.regime=trop.simmap.allbirds.i,trim.class="fruit")
		regime.class.df<-return.class.df(out$regime.simmap,out$regime.class.object)

		#set lineages not in "fruit" subgroup to evolve via BM
		regime.class.df[,which(colnames(out$regime.simmap$mapped.edge)=='Z')+1]=1

		#trim the fruit region simmap to just lineages that are fruit at the tips, adjust class.df accordingly
		regime.simmap.region.trimmed<-drop.tip.simmap(out$regime.simmap,out$regime.simmap$tip.label[which(!out$regime.simmap$tip.label%in%names(M))])
		regime.class.df.trimmed<-regime.class.df[,c(1,match(colnames(regime.simmap.region.trimmed$mapped.edge),colnames(out$regime.simmap$mapped.edge))+1)]		

		regime.simmap.region.root=max(nodeHeights(out$regime.simmap))
		regime.simmap.region.trimmed.root=max(nodeHeights(regime.simmap.region.trimmed))

		###	adjust the class.df and class.object if the trimmed tree has a younger root than the ancestral tree
		if(round(regime.simmap.region.root,5)!=round(regime.simmap.region.trimmed.root,5)){

			trimmed.class.object<-CreateClassObject(regime.simmap.region.trimmed)
			shifted.times<-trimmed.class.object$times+(regime.simmap.region.root-regime.simmap.region.trimmed.root)
			new.regime.class.df.trimmed<-regime.class.df.trimmed[c(which(round(out$regime.class.object$times,5)==round(min(shifted.times),5)):dim(regime.class.df.trimmed)[1]),]
			
			new.regime.class.df.trimmed$interval<-c(1:dim(new.regime.class.df.trimmed)[1])
			out$regime.class.object$times<-out$regime.class.object$times[which(round(out$regime.class.object$times,5)>=round(min(shifted.times),5))]-round(regime.simmap.region.root-regime.simmap.region.trimmed.root,5)
			#forces time to start at root of trimmed tree; would be better to pass times directly to new_list_function to avoid overwriting this slot of the out object, which could lead to errors
		}
		
		new_list_function<-create.function.list(regime.simmap.region.trimmed,out$regime.class.object,regime.class.df.trimmed)
		o3B<-fit_t_general(regime.simmap.region.trimmed, M, fun=new_list_function,diagnostic=T,echo=T,constraint=FALSE)			 
		DDexpgeo_mvMORPH.lnL<-o3B$LogLik
		DDexpgeo_mvMORPH.sig2<-o3B$rates[2,1]
		DDexpgeo_mvMORPH.r<-o3B$rates[1,1]
		DDexpgeo_mvMORPH.z0<-o3B$anc
		DDexpgeo_mvMORPH.conv<-o3B$convergence

##STOPPED here		
		
		int<-c(colnames(master)[90+j],"fruit",length(M),o1$LogLik,o1$AIC,o1$sigma,o1$theta,o1$convergence,o2$LogLik,o2$AIC,o2$sigma,o2$alpha,o2$theta,o2$convergence,iter,o3$LogLik,o3$rates[2,1],o3$rates[1,1],o3$anc,o3$convergence)
		print(int)
		res.mat[j,]<-int
		write.csv(res.mat,file="modelfits_frugivores_firstmap.csv")
		print(j)
	}	

write.csv(res.mat,file="modelfits_frugivores_firstmap.csv")

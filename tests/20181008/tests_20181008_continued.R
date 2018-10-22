##This script simulates multi-slope DD model to test with subgroup pruning

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
source('CreateClassbyClassObject_mvMORPH.R')
source('CreateSMatrix_subgroup.R')
source('sim_DDmulti_3S.R')

#source('~/ddexp/R/fit_t_general.R')
#source('~/ddexp/R/generalized_functions.R')
#source('~/ddexp/R/CreateGeobyClassObject_mvMORPH.R')
#source('~/ddexp/R/make.simmap.BGB.R')
#source('~/Dropbox/Scripts/R scripts/trimSimmap.R')
#source('~/Dropbox/Scripts/R scripts/CreateClassObject.R')
#source('~/ddexp/R/CreateClassbyClassObject_mvMORPH.R')
#source('~/ddexp/tests/20181008/CreateSMatrix_subgroup.R', chdir = TRUE)
#source('~/ddexp/tests/20181008/sim_DDmulti_3S.R', chdir = TRUE)
#
master<-read.csv("BIRDS_MASTER.csv")
tree<-read.tree("treefile_allbirds_simple.tree")

load('diet.simmap.allbirds.RData')
load('all.simmap.trop.RData')

#master<-read.csv("~/Dropbox/Temperate:Tropical/BIRDS_MASTER.csv")
#tree<-read.tree("~/Dropbox/Temperate:Tropical/biogeography_data/Non-stratified simple DEC allbirds/treefile_allbirds_simple.tree")
#
#load('~/Dropbox/Temperate:Tropical/simmaps/diet.simmap.allbirds.RData')
#load('~/Dropbox/Temperate:Tropical/simmaps/all.simmap.trop.RData')


#diet.simmaps were generated with all data, whereas biogeography was limited to a subset of species; need to trim to reflect this
diet.simmap.allbirds.i<-diet.simmap.allbirds[[3]]
diet.simmap.allbirds.i<-drop.tip.simmap(diet.simmap.allbirds.i,diet.simmap.allbirds.i$tip.label[which(!diet.simmap.allbirds.i$tip.label%in%tree$tip.label)])

trop.simmap.allbirds.i<-all.simmap.trop[[4]]
trop.simmap.allbirds.i<-drop.tip.simmap(trop.simmap.allbirds.i,trop.simmap.allbirds.i$tip.label[which(!trop.simmap.allbirds.i$tip.label%in%tree$tip.label)])

nectar<-subset(master,diet=="nectar")

region<-subset(master,H==1)
region.tree<-drop.tip(tree,tree$tip.label[which(!tree$tip.label%in%region$Species_name_from_birdtree)])
nectar<-nectar[nectar$Species_name_from_birdtree%in%diet.simmap.allbirds.i$tip.label,]

j=1
M<-nectar[,90+j]
names(M)<-nectar[,53]
M<-subset(M,M!="NA" & names(M) %in% region.tree$tip.label)
nc<-name.check(region.tree,M)
#nc
if(is.list(nc)){
	subtree<-drop.tip(region.tree,nc$tree_not_data)
	} else{
	subtree<-region.tree
	}
M<-M[subtree$tip.label]

### NOW fitting models that vary according to tropical/temperate membership

#trim diet simmap to species found in focal region
diet.simmap.region<-drop.tip.simmap(diet.simmap.allbirds.i,diet.simmap.allbirds.i$tip.label[which(!diet.simmap.allbirds.i$tip.label%in%region.tree$tip.label)])

trop.simmap.trimmed<-drop.tip.simmap(trop.simmap.allbirds.i,trop.simmap.allbirds.i$tip.label[which(!trop.simmap.allbirds.i$tip.label%in%subtree$tip.label)])

##prepare simmap and fit DDexp to subgroup *regime

#need to build a class.df where diversity represents the number of species in a REGION that are ALSO in the nectar category

out<-try(CreateClassbyClassObject_mvMORPH(map.guild=diet.simmap.region,map.regime=trop.simmap.allbirds.i,trim.class="nectar"))
if(class(out)=="try-error"){out<-CreateClassbyClassObject_mvMORPH(map.guild=diet.simmap.region,map.regime=trop.simmap.allbirds.i,trim.class="nectar",rnd=6)}

slope.list=list("1","5","Z")
smat<-CreateSMatrix_subgroup(out$regime.class.object,S.cats=slope.list)


##conduct simulation study:

pars<-rbind(expand.grid(rep(0.1,50),-0.05,c(0.05,0,-0.025,-0.05,-.1),0,0),expand.grid(rep(0.1,50),c(0.05,0,-0.025,-0.05,-.1),-0.05,0,0))
pars[,1]<-0.01/(0.5*(exp(pars[,2]*86)+exp(pars[,3]*33))) #keep sig2 at tips constant

#simdata<-sim_multiDD_subgroup(out$regime.simmap,pars=pars, S.matrix=smat,rnd=5,min.Nsegments=50000)
#save(simdata,file="simdata.RData")
load('simdata.RData')
##then write script to fit everything

regime.class.df<-return.class.df(out$regime.simmap,out$regime.class.object)

#set lineages not in "nectar" subgroup to evolve via BM
regime.class.df[,which(colnames(out$regime.simmap$mapped.edge)=='Z')+1]=1

#trim the nectar region simmap to just lineages that are nectar at the tips, adjust class.df accordingly
regime.simmap.region.trimmed<-drop.tip.simmap(out$regime.simmap,out$regime.simmap$tip.label[which(!out$regime.simmap$tip.label%in%names(M))])


regime.class.df.trimmed<-regime.class.df[,c(1,match(colnames(regime.simmap.region.trimmed$mapped.edge),colnames(out$regime.simmap$mapped.edge))+1)]		

regime.simmap.region.root=max(nodeHeights(out$regime.simmap))
regime.simmap.region.trimmed.root=max(nodeHeights(regime.simmap.region.trimmed))

if(round(regime.simmap.region.root,5)!=round(regime.simmap.region.trimmed.root,5)){

	trimmed.class.object<-CreateClassObject(regime.simmap.region.trimmed,rnd=5)
	shifted.times<-trimmed.class.object$times+(regime.simmap.region.root-regime.simmap.region.trimmed.root)
	new.regime.class.df.trimmed<-regime.class.df.trimmed[c(which(round(out$regime.class.object$times,5)==round(min(shifted.times),5)):dim(regime.class.df.trimmed)[1]),]
	
	new.regime.class.df.trimmed$interval<-c(1:dim(new.regime.class.df.trimmed)[1])
	out$regime.class.object$times<-out$regime.class.object$times[which(round(out$regime.class.object$times,5)>=round(min(shifted.times),5))]-round(regime.simmap.region.root-regime.simmap.region.trimmed.root,5)
	#forces time to start at root of trimmed tree; would be better to pass times directly to new_list_function to avoid overwriting this slot of the out object, which could lead to errors
	regime.class.df.trimmed<-new.regime.class.df.trimmed
}

new_list_function<-create.function.list(regime.simmap.region.trimmed,out$regime.class.object,regime.class.df.trimmed)



res.mat<-matrix(ncol=68,nrow=dim(simdata)[1])
colnames(res.mat)<-c('sim.sig2_1','sim.sig2_2','sim.sig2_3','sim.r1','sim.r2','sim.r3','sim.z0','est.sig2_1','est.sig2_2','est.sig2_3','est.r1','est.r2','est.r3','est.z0',"BM1.lnL","BM1.AICc","BM1.sig","BM1.z0","BM1.convergence","BM1.hess.value","OU1.lnL","OU1.AICc","OU1.sig","OU1.alpha","OU1.z0","OU1.convergence","OU1.hess.value","DDexp.lnL","DDexp.sig","DDexp.r","DDexp.z0","DDexp.convergence","DDexp.hess.value","BMM.lnL","BMM.AICc","BMM.cat_1","BMM.cat_2","BMM.sig_1","BMM.sig_2","BMM.z0","BMM.convergence","BMM.hess.value","OUM.lnL","OUM.AICc","OUM.sig","OUM.alpha","OUM.cat_1","OUM.cat_2","OUM.z0_1","OUM.z0_2","OUM.convergence","OUM.hess.value","DDM.lnL","DDM.AICc","DDM.convergence","DDM.hess.value","delAICc.BM1","delAICc.OU1","delAICc.DD1","delAICc.BMM","delAICc.OUM","delAICc.DDM","BM1.wi","OU1.wi","DD1.wi","BMM.wi","OUM.wi","DDM.wi")


for(i in 138:dim(simdata)[1]){
		M<-simdata[i,]
		names(M)<-colnames(simdata)
		
		#pull out species in subgroup only
		M<-M[which(names(M)%in%subtree$tip.label)]
		
		o1<-try(mvBM(subtree,M,model=c("BM1"),optimization="subplex"))
		if(class(o1)[1]=="try-error"){
			o1<-mvBM(subtree,M,model=c("BM1"))
			if(o1$convergence!=0 | o1$hess.values!=0){
			o1<-mvBM(subtree,M,model=c("BM1"),optimization="Nelder-Mead")
			}
		} else {
			if(o1$convergence!=0 | o1$hess.values!=0){
			o1<-mvBM(subtree,M,model=c("BM1"),optimization="Nelder-Mead")
			}
		}
	
		o2<-try(mvOU(subtree,M,model=c("OU1"),optimization="subplex"))
		if(class(o2)[1]=="try-error"){
			o2<-try(mvOU(subtree,M,model=c("OU1")))
			if((class(o2)[1]=="try-error") || (o2$convergence!=0 | o2$hess.values!=0)){
			o2<-mvOU(subtree,M,model=c("OU1"),optimization="Nelder-Mead")
			}
		} else {
			if(o2$convergence!=0 | o2$hess.values!=0){
			o2<-try(mvOU(subtree,M,model=c("OU1"),optimization="Nelder-Mead"))
			if(class(o2)[1]=="try-error"){
				o2<-mvOU(subtree,M,model=c("OU1"))
				}
			}
		}

		##prepare simmap and fit DDexp to subgroup
				
		#trim this diet*region simmap to contain only branches in "nectar" subgroup (see Drury et al. 2018 PLOS Biol for description of algorithm)
		nectar.simmap.region<-trimSimmap(diet.simmap.region,trim.class="nectar")
		
		#create 'class object' for subsequent scripts
		out<-try(CreateClassObject(nectar.simmap.region))
		if(class(out)=="try-error"){out<-CreateClassObject(nectar.simmap.region,rnd=6)}
		
		#create class.df that has the diversity of lineages in each simmap class
		nectar.class.df<-return.class.df_subgroup(nectar.simmap.region,out)
		
		#set lineages not in "nectar" subgroup to evolve via BM
		nectar.class.df[,which(colnames(nectar.simmap.region$mapped.edge)!='nectar')+1]=1 
	
		#trim the nectar region simmap to just lineages that are nectar at the tips, adjust class.df accordingly
		nectar.simmap.region.trimmed<-drop.tip.simmap(nectar.simmap.region,nectar.simmap.region$tip.label[which(!nectar.simmap.region$tip.label%in%names(M))])

		nectar.class.df.trimmed<-nectar.class.df[,c(1,match(colnames(nectar.simmap.region.trimmed$mapped.edge),colnames(nectar.simmap.region$mapped.edge))+1)]		
	
		nectar.simmap.region.root=max(nodeHeights(nectar.simmap.region))
		nectar.simmap.region.trimmed.root=max(nodeHeights(nectar.simmap.region.trimmed))

		###	adjust the class.df and class.object if the trimmed tree has a younger root than the ancestral tree
		if(round(nectar.simmap.region.root,5)!=round(nectar.simmap.region.trimmed.root,5)){
			
			trimmed.class.object<-CreateClassObject(nectar.simmap.region.trimmed)
			shifted.times<-trimmed.class.object$times+(nectar.simmap.region.root-nectar.simmap.region.trimmed.root)
			
			new.nectar.class.df.trimmed<-nectar.class.df.trimmed[c(which(round(out$times,5)==round(min(shifted.times),5)):dim(nectar.class.df.trimmed)[1]),]
			
			new.nectar.class.df.trimmed$interval<-c(1:dim(new.nectar.class.df.trimmed)[1])
			
			nectar.class.df.trimmed<-new.nectar.class.df.trimmed
			out$times<-out$times[which(round(out$times,5)>=round(min(shifted.times),5))]-round(nectar.simmap.region.root-nectar.simmap.region.trimmed.root,5)
			#forces time to start at root of trimmed tree; would be better to pass times directly to new_list_function to avoid overwriting this slot of the out object, which could lead to errors
			regime.class.df.trimmed<-new.nectar.class.df.trimmed
		}
		
		list_function<-create.function.list(nectar.simmap.region.trimmed,out,nectar.class.df.trimmed)
		o3<-fit_t_general(nectar.simmap.region.trimmed, M, fun=list_function,diagnostic=T,echo=T)			 
		if(o3$convergence!=0 | o3$hess.values!=0){
		o3<-fit_t_general(nectar.simmap.region.trimmed, M, fun=list_function,diagnostic=T,echo=T,method="BB")			 
		}
		
		o1B<-try(mvBM(trop.simmap.trimmed,M,model=c("BMM"),optimization="subplex"))
		if(class(o1B)[1]=="try-error"){
			o1B<-mvBM(trop.simmap.trimmed,M,model=c("BMM"))
			if(o1B$convergence!=0 | o1B$hess.values!=0){
			o1B<-mvBM(trop.simmap.trimmed,M,model=c("BMM"),optimization="Nelder-Mead")
			}
		} else {
			if(o1B$convergence!=0 | o1B$hess.values!=0){
			o1B<-mvBM(trop.simmap.trimmed,M,model=c("BMM"),optimization="Nelder-Mead")
			}
		}
		o2B<-try(mvOU(trop.simmap.trimmed,M,model=c("OUM"),optimization="subplex"))
		if(class(o2B)=="try-error"){
			o2B<-mvOU(trop.simmap.trimmed,M,model=c("OUM"))
			if(o2B$convergence!=0 | o2B$hess.values!=0){
			o2B<-mvOU(trop.simmap.trimmed,M,model=c("OUM"),optimization="Nelder-Mead")
			}
		} else {
			if(o2B$convergence!=0 | o2B$hess.values!=0){
			o2B<-try(mvOU(trop.simmap.trimmed,M,model=c("OUM"),optimization="Nelder-Mead"))
			if(class(o2B)=="try-error"){
			o2B<-mvOU(trop.simmap.trimmed,M,model=c("OUM"))
			}
			}
		}

		o3B<-fit_t_general(regime.simmap.region.trimmed, M, fun=new_list_function,diagnostic=T,echo=T,constraint=FALSE)	
		if(o3B$convergence!=0 | o3B$hess.values!=0){
		o3B<-fit_t_general(regime.simmap.region.trimmed, M, fun=new_list_function,diagnostic=T,echo=T,constraint=FALSE,method="BB",control=list(maxit=50000,maxfeval=50000))				 		
		}
		
		delAICc.BM1<-o1$AICc-min(o1$AICc,o2$AICc,o3$AICc,o1B$AICc,o2B$AICc,o3B$AICc)
		delAICc.OU1<-o2$AICc-min(o1$AICc,o2$AICc,o3$AICc,o1B$AICc,o2B$AICc,o3B$AICc) 
		delAICc.DD1<-o3$AICc-min(o1$AICc,o2$AICc,o3$AICc,o1B$AICc,o2B$AICc,o3B$AICc)
		delAICc.BMM<-o1B$AICc-min(o1$AICc,o2$AICc,o3$AICc,o1B$AICc,o2B$AICc,o3B$AICc) 
		delAICc.OUM<-o2B$AICc-min(o1$AICc,o2$AICc,o3$AICc,o1B$AICc,o2B$AICc,o3B$AICc)
		delAICc.DDM<-o3B$AICc-min(o1$AICc,o2$AICc,o3$AICc,o1B$AICc,o2B$AICc,o3B$AICc) 
		
		all=sum(exp(-0.5*delAICc.BM1),exp(-0.5*delAICc.OU1),exp(-0.5*delAICc.DD1),exp(-0.5*delAICc.BMM),exp(-0.5*delAICc.OUM),exp(-0.5*delAICc.DDM))
		BM1.wi<-exp(-0.5*delAICc.BM1)/all
		OU1.wi<-exp(-0.5*delAICc.OU1)/all
		DD1.wi<-exp(-0.5*delAICc.DD1)/all
		BMM.wi<-exp(-0.5*delAICc.BMM)/all
		OUM.wi<-exp(-0.5*delAICc.OUM)/all
		DDM.wi<-exp(-0.5*delAICc.DDM)/all

		
		int<-c(as.numeric(rep(pars[i,1],3)),as.numeric(pars[i,2:5]),o3B$rates[,which(colnames(o3B$rates)==slope.list[[1]])][2],o3B$rates[,which(colnames(o3B$rates)==slope.list[[2]])][2],o3B$rates[,which(colnames(o3B$rates)==slope.list[[3]])][2],o3B$rates[,which(colnames(o3B$rates)==slope.list[[1]])][1],o3B$rates[,which(colnames(o3B$rates)==slope.list[[2]])][1],o3B$rates[,which(colnames(o3B$rates)==slope.list[[3]])][1],o3B$anc,o1$LogLik,o1$AICc,o1$sigma,o1$theta,o1$convergence,o1$hess.value,o2$LogLik,o2$AICc,o2$sigma,o2$alpha,o2$theta,o2$convergence,o2$hess.value,o3$LogLik,o3$rates[2,1],o3$rates[1,1],o3$anc,o3$convergence,o3$hess.value,o1B$LogLik,o1B$AICc,dimnames(o1B$sigma)[[3]][1],dimnames(o1B$sigma)[[3]][2],o1B$sigma[1],o1B$sigma[2],o1B$theta,o1B$convergence,o1B$hess.value,o2B$LogLik,o2B$AICc,o2B$sigma,o2B$alpha,rownames(o2B$theta)[1],rownames(o2B$theta)[2],o2B$theta[1],o2B$theta[2],o2B$convergence,o2B$hess.value,o3B$LogLik,o3B$AICc,o3B$convergence,o3B$hess.value,delAICc.BM1,delAICc.OU1,delAICc.DD1,delAICc.BMM,delAICc.OUM,delAICc.DDM,BM1.wi,OU1.wi,DD1.wi,BMM.wi,OUM.wi,DDM.wi)
		res.mat[i,]<-int
		print(int)
		print(i)
		write.csv(res.mat,file="ddm_simulation_study_20181012.csv")
	}	

write.csv(res.mat,file="ddm_simulation_study_20181012.csv")


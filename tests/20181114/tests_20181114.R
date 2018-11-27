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

source('fit_t_general_options.R')
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



res.mat<-matrix(ncol=24,nrow=dim(simdata)[1])
colnames(res.mat)<-c('sim.sig2_1','sim.sig2_2','sim.sig2_3','sim.r1','sim.r2','sim.r3','sim.z0','est.sig2_1','est.sig2_2','est.sig2_3','est.r1','est.r2','est.r3','est.z0',"BM1.lnL","BM1.AICc","BM1.sig","BM1.z0","BM1.convergence","BM1.hess.value","DDM_constraint.lnL","DDM_constraint.AICc","DDM_constraint.convergence","DDM_constraint.hess.value")


for(i in 1:dim(simdata)[1]){
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
			
		beta.constraint<-vector(length=3)
		sigma.constraint<-c(1,1,1)
		beta.constraint[which(colnames(regime.simmap.region.trimmed$mapped.edge)=="Z")]<-NA
		beta.constraint[which(colnames(regime.simmap.region.trimmed$mapped.edge)!="Z")]<-1:2
		
		o3B<-fit_t_general(regime.simmap.region.trimmed, M, fun=new_list_function,diagnostic=T,echo=T,constraint=list(sigma=sigma.constraint, beta=beta.constraint))	
		if(o3B$convergence!=0 | o3B$hess.values!=0){
		o3B<-fit_t_general(regime.simmap.region.trimmed, M, fun=new_list_function,diagnostic=T,echo=T,constraint=list(sigma=sigma.constraint, beta=beta.constraint),method="BB",control=list(maxit=50000,maxfeval=50000))				 		
		}
		
		int<-c(as.numeric(rep(pars[i,1],3)),as.numeric(pars[i,2:5]),o3B$rates[,which(colnames(o3B$rates)==slope.list[[1]])][2],o3B$rates[,which(colnames(o3B$rates)==slope.list[[2]])][2],o3B$rates[,which(colnames(o3B$rates)==slope.list[[3]])][2],o3B$rates[,which(colnames(o3B$rates)==slope.list[[1]])][1],o3B$rates[,which(colnames(o3B$rates)==slope.list[[2]])][1],o3B$rates[,which(colnames(o3B$rates)==slope.list[[3]])][1],o3B$anc,o1$LogLik,o1$AICc,o1$sigma,o1$theta,o1$convergence,o1$hess.value,o3B$LogLik,o3B$AICc,o3B$convergence,o3B$hess.value)
		res.mat[i,]<-int
		print(int)
		print(i)
		write.csv(res.mat,file="ddm_simulation_study_20181114.csv")
	}	

write.csv(res.mat,file="ddm_simulation_study_20181114.csv")


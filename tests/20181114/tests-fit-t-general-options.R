# Tests (from Julien Clavel)
source('~/ddexp/R/fit_t_general_options.R')
library(mvMORPH)
library(RPANDA)

# Simulated dataset
set.seed(14)
# Generating a random tree
tree<-pbtree(n=50)

# Setting the regime states of tip species
sta<-as.vector(c(rep("Forest",20),rep("Savannah",30))); names(sta)<-tree$tip.label

# Making the simmap tree with mapped states
tree<-make.simmap(tree,sta , model="ER", nsim=1)
col<-c("blue","orange"); names(col)<-c("Forest","Savannah")

# Plot of the phylogeny for illustration
plotSimmap(tree,col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)

# trait
trait <- rTraitCont(tree)

# function
data(InfTemp)
require(pspline)
spline_result <- sm.spline(x=InfTemp[,1],y=InfTemp[,2], df=50)
env_func <- function(t){predict(spline_result,t)}
t<-unique(InfTemp[,1])

# We build the interpolated smoothing spline function
env_data<-splinefun(t,env_func(t))
funEnv <- list(env_data, env_data)

# tests
# The unconstrained model
fit_t_general(tree, trait, fun=funEnv, error=NULL, beta=NULL, sigma=NULL, 
              model=c("exponential","linear"), method=c("L-BFGS-B","BB"), 
              upper=Inf, lower=-20, control=list(maxit=20000), diagnostic=TRUE, 
              echo=TRUE, constraint=NULL)

# constrain only sigma
fit_t_general(tree, trait, fun=funEnv, error=NULL, beta=NULL, sigma=NULL, 
              model=c("exponential","linear"), method=c("L-BFGS-B","BB"), 
              upper=Inf, lower=-20, control=list(maxit=20000), diagnostic=TRUE, 
              echo=TRUE, constraint=list(sigma=c(1,1)))

# constrain one of the beta to be zero and common sigma
fit_t_general(tree, trait, fun=funEnv, error=NULL, beta=NULL, sigma=NULL, 
              model=c("exponential","linear"), method=c("L-BFGS-B","BB"), 
              upper=Inf, lower=-20, control=list(maxit=20000), diagnostic=TRUE, 
              echo=TRUE, constraint=list(sigma=c(1,1), beta=c(NA,1)))

# Equivalent to the constraint model from the previous version:
fit_t_general(tree, trait, fun=funEnv, error=NULL, beta=NULL, sigma=NULL, 
              model=c("exponential","linear"), method=c("L-BFGS-B","BB"), 
              upper=Inf, lower=-20, control=list(maxit=20000), diagnostic=TRUE, 
              echo=TRUE, constraint=list(sigma=c(1,1), beta=c(1,1)))


###ADDING IN a test using a simple case

require(deSolve)
require(geiger)
require(phytools)
library(corpcor)
require(methods)
require(tools)
require(RPANDA)
require(data.table)
require(mvMORPH)

source('~/ddexp/R/generalized_functions.R')
source('~/ddexp/R/CreateGeobyClassObject_mvMORPH.R')
source('~/ddexp/R/make.simmap.BGB.R')
source('~/Dropbox/Scripts/R scripts/trimSimmap.R')
source('~/Dropbox/Scripts/R scripts/CreateClassObject.R')
source('~/ddexp/R/CreateClassbyClassObject_mvMORPH.R')

master<-read.csv("~/Dropbox/Temperate:Tropical/BIRDS_MASTER.csv")
tree<-read.tree("~/Dropbox/Temperate:Tropical/biogeography_data/Non-stratified simple DEC allbirds/treefile_allbirds_simple.tree")

load('~/Dropbox/Temperate:Tropical/simmaps/diet.simmap.allbirds.RData')
load('~/Dropbox/Temperate:Tropical/simmaps/all.simmap.trop.RData')

#diet.simmaps were generated with all data, whereas biogeography was limited to a subset of species; need to trim to reflect this
diet.simmap.allbirds.i<-diet.simmap.allbirds[[1]]
diet.simmap.allbirds.i<-drop.tip.simmap(diet.simmap.allbirds.i,diet.simmap.allbirds.i$tip.label[which(!diet.simmap.allbirds.i$tip.label%in%tree$tip.label)])

trop.simmap.allbirds.i<-all.simmap.trop[[1]]
trop.simmap.allbirds.i<-drop.tip.simmap(trop.simmap.allbirds.i,trop.simmap.allbirds.i$tip.label[which(!trop.simmap.allbirds.i$tip.label%in%tree$tip.label)])

region.vec<-c("AB","CD","E","F","G","H","I","J","L")

res.mat<-matrix(nrow=7*length(region.vec),ncol=69)
colnames(res.mat)<-c("region","trait","diet.category","no","no.lineages.in.each.regime","BM1.lnL","BM1.AICc","BM1.sig","BM1.z0","BM1.convergence","BM1.hess.value","OU1.lnL","OU1.AICc","OU1.sig","OU1.alpha","OU1.z0","OU1.convergence","OU1.hess.value","DDexp.lnL","DDexp.sig","DDexp.r","DDexp.z0","DDexp.convergence","DDexp.hess.value","BMM.lnL","BMM.AICc","BMM.cat_1","BMM.cat_2","BMM.sig_1","BMM.sig_2","BMM.z0","BMM.convergence","BMM.hess.value","OUM.lnL","OUM.AICc","OUM.sig","OUM.alpha","OUM.cat_1","OUM.cat_2","OUM.z0_1","OUM.z0_2","OUM.convergence","OUM.hess.value","DDM.lnL","DDM.AICc","DDM.cat_1","DDM.cat_2","DDM.cat_3","DDM.sig_1","DDM.sig_2","DDM.sig_3","DDM.rate_1","DDM.rate_2","DDM.rate_3","DDM.z0","DDM.convergence","DDM.hess.value","delAICc.BM1","delAICc.OU1","delAICc.DD1","delAICc.BMM","delAICc.OUM","delAICc.DDM","BM1.wi","OU1.wi","DD1.wi","BMM.wi","OUM.wi","DDM.wi")

fruit<-subset(master,diet=="fruit")

i=1
j=2
	
#trim tree down to region of interest

if(i==1){
region<-subset(master,A==1 | B==1)
region.tree<-drop.tip(tree,tree$tip.label[which(!tree$tip.label%in%region$Species_name_from_birdtree)])
}
##need to further trim until simmap is constructed for expanded set of species
fruit<-fruit[fruit$Species_name_from_birdtree%in%diet.simmap.allbirds.i$tip.label,]

M<-fruit[,90+j]
names(M)<-fruit[,53]
M<-subset(M,M!="NA" & names(M) %in% region.tree$tip.label)
nc<-name.check(region.tree,M)
#nc
if(is.list(nc)){
	subtree<-drop.tip(region.tree,nc$tree_not_data)
	} else{
	subtree<-region.tree
	}
M<-M[subtree$tip.label]

#trim diet simmap to species found in focal region
diet.simmap.region<-drop.tip.simmap(diet.simmap.allbirds.i,diet.simmap.allbirds.i$tip.label[which(!diet.simmap.allbirds.i$tip.label%in%region.tree$tip.label)])

#need to build a class.df where diversity represents the number of species in a REGION that are ALSO in the fruit category

out<-try(CreateClassbyClassObject_mvMORPH(map.guild=diet.simmap.region,map.regime=trop.simmap.allbirds.i,trim.class="fruit"))
if(class(out)=="try-error"){out<-CreateClassbyClassObject_mvMORPH(map.guild=diet.simmap.region,map.regime=trop.simmap.allbirds.i,trim.class="fruit",rnd=6)}
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

	trimmed.class.object<-CreateClassObject(regime.simmap.region.trimmed,rnd=5)
	shifted.times<-trimmed.class.object$times+(regime.simmap.region.root-regime.simmap.region.trimmed.root)
	new.regime.class.df.trimmed<-regime.class.df.trimmed[c(which(round(out$regime.class.object$times,5)==round(min(shifted.times),5)):dim(regime.class.df.trimmed)[1]),]
	
	new.regime.class.df.trimmed$interval<-c(1:dim(new.regime.class.df.trimmed)[1])
	out$regime.class.object$times<-out$regime.class.object$times[which(round(out$regime.class.object$times,5)>=round(min(shifted.times),5))]-round(regime.simmap.region.root-regime.simmap.region.trimmed.root,5)
	#forces time to start at root of trimmed tree; would be better to pass times directly to new_list_function to avoid overwriting this slot of the out object, which could lead to errors
	regime.class.df.trimmed<-new.regime.class.df.trimmed
}

new_list_function<-create.function.list(regime.simmap.region.trimmed,out$regime.class.object,regime.class.df.trimmed)

#as before
o3B<-fit_t_general(regime.simmap.region.trimmed, M, fun=new_list_function,diagnostic=T,echo=T,constraint=NULL)
	
#with desired constraint
o3B<-fit_t_general(regime.simmap.region.trimmed, M, fun=new_list_function,diagnostic=T,echo=T,constraint=list(sigma=c(1,1,1), beta=c(NA,1,2)))	


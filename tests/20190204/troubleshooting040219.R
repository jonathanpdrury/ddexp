##testing two slope version

set.seed(1234)

fiftytipyuletrees<-pbtree(n=50,nsim=1)
regions<-c("trop","trop","trop","temp","temp")
phylo<-fiftytipyuletrees
samp.reg<-sample(regions)
trait<-c(rep(samp.reg[1],10),rep(samp.reg[2],10),rep(samp.reg[3],10),rep(samp.reg[4],10),rep(samp.reg[5],10))
names(trait)<-phylo$tip.label
smap<-make.simmap(phylo,trait,model="ARD")
source('~/Dropbox/Temperate:Tropical/model_development/DD multi-slope models/sim_DDmulti.R')



source('~/Dropbox/Temperate:Tropical/model_development/DD multi-slope models/DDexpMulti_nogeo_ADiag.R')
source('~/Dropbox/Temperate:Tropical/model_development/DD multi-slope models/DDlinMulti_nogeo_ADiag.R')

source('~/Dropbox/Scripts/R scripts/resortSMatrix.R')
source('~/Dropbox/Scripts/R scripts/CreateSMatrix.R')

r.matrix= CreateSMatrix(CreateClassObject(smap),S.cats=c("trop","temp"))
trait<-try(sim_multiDD(smap,pars=matrix(c(0.05,-0.05,-0.1,0),nrow=1),plot=F,S.matrix=r.matrix,rnd=5,model="exponential"))

rmats<-resortSMatrix(phylo,r.matrix)
ddm.ob<-createModel_DDexp_multi(phylo,rmats)
params0<-c(0,log(sqrt(var(trait[1,])/max(nodeHeights(phylo)))),-1e-4,-1e-4)	
o1a<-fitTipData(ddm.ob,trait[1,],params0,GLSstyle=TRUE)

#$value
#[1] -36.6613
#
#$inferredParams
#         m0   logsigma0          r1          r2 
#-0.24528295 -1.24126464 -0.07943798 -0.10991249 
#
#$convergence
#[1] 0

ddl.ob<-createModel_DDlin_multi(phylo,rmats)
o1al<-fitTipData(ddl.ob,trait[1,],params0,GLSstyle=TRUE)

#> o1al
#$value
#[1] -34.60216
#
#$inferredParams
#          m0    logsigma0           r1           r2 
#-0.245224142 -1.577608135 -0.001132290 -0.001420973 
#
#$convergence
#[1] 0

source('~/ddexp/R/fit_t_DD_BETA.R')


o1b<-fit_t_DD(phylo,data=trait[1,],model="exponential",regime.map=smap,method=method)
#> o1b
#$LogLik
#[1] 36.67312
#
#$AIC
#[1] -65.34624
#
#$AICc
#[1] -64.45735
#
#$rates
#             temp        trop
#beta  -0.11002697 -0.07952908
#sigma  0.08367657  0.08367657
#
#$anc
#[1] -0.2458498
#
#$convergence
#[1] 0
#
#$hess.values
#[1] 0
#
#$error
#NULL
#
#$param
#[1] -0.11002697 -0.07952908 -2.48079623
#
#$phyloTrans
#NULL
#
#> exp(-1.24126464)^2
#[1] 0.08353168

-getDataLikelihood(ddm.ob, trait[1,], params=c(o1b$anc,log(sqrt(o1b$rates[2,1])),o1b$rates[1,2],o1b$rates[1,1]))
#[1] 36.66129




o1bl<-fit_t_DD(phylo,data=trait[1,],model="linear",regime.map=smap,method="BB")
#$LogLik
#[1] 31.72709
#
#$AIC
#[1] -55.45417
#
#$AICc
#[1] -54.56528
#
#$rates
#               temp          trop
#beta  -0.0004840021 -0.0004616578
#sigma  0.0242001053  0.0242001053
#
#$anc
#[1] -0.2450647
#
#$convergence
#[1] 0
#
#$hess.values
#[1] 0
#
#$error
#NULL
#
#$param
#[1] -0.0004840021 -0.0004616578  0.0242001053
#
#$phyloTrans
#NULL

-getDataLikelihood(ddl.ob, trait[1,], params=c(o1bl$anc,log(sqrt(o1bl$rates[2,1])),o1bl$rates[1,2],o1bl$rates[1,1]))
#[1] 31.72649

##these are returning the same likelihood, but for whatever reason fit_t_DD is not optimising correctly


###two-slope version is working fine for exponential case, but not linear case, where optimisation doesn't reach maximum

######################################
######################################
######################################
######################################

###now back to one slope version

o1a<-fit_t_comp(phylo,trait[1,],model="DDexp")
o1b<-fit_t_DD(phylo,data=trait[1,],model="exponential",method=method)

source('~/Dropbox/Scripts/R scripts/DDexp_nogeo_ADiag.R')
ddm.ob<-createModel_DDexp(phylo)
params0<-c(0,log(sqrt(var(trait[1,])/max(nodeHeights(phylo)))),-1e-4)	
o1a<-fitTipData(ddm.ob,trait[1,],params0,GLSstyle=TRUE)
-getDataLikelihood(ddm.ob, trait[1,], params=c(o1b$anc,log(sqrt(o1b$rates[2,1])),o1b$rates[1,1]))

##these are returning the same likelihood, but for whatever reason fit_t_DD is not optimising correctly for either model

######################################
######################################
######################################
######################################

#one slope with biogeography

library(RPANDA)

data(BGB.examples)
##Example with a non-stratified tree
Canidae.geography.object<-CreateGeoObject_BioGeoBEARS(full.phylo=BGB.examples$Canidae.phylo,
ana.events=BGB.examples$Canidae.ana.events, clado.events=BGB.examples$Canidae.clado.events)
data<-rnorm(length(BGB.examples$Canidae.phylo$tip.label))
names(data)<-BGB.examples$Canidae.phylo$tip.label

DDlin.geo.fit<-fit_t_comp(BGB.examples$Canidae.phylo, data,model="DDlin", geography.object=Canidae.geography.object)
DDexp.geo.fit<-fit_t_comp(BGB.examples$Canidae.phylo, data, model="DDexp", geography.object=Canidae.geography.object)


###checking fit_t_DD
source('~/ddexp/R/make.simmap.BGB.R', chdir = TRUE)

require(geiger)
require(data.table)
geo.map=make.simmap.BGB(BGB.examples$Canidae.phylo,BGB.examples$Canidae.phylo,ana.events=BGB.examples$Canidae.ana.events,clado.events=BGB.examples$Canidae.clado.events,rnd=6)$geo.simmap
fit_t_DD(BGB.examples$Canidae.phylo,data,model="exponential",geo.map=geo.map)

######################################
######################################
######################################
######################################

#one slope with subgroup pruning (see tests_20180928.R)

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
source('~/Dropbox/Temperate:Tropical/model_development/DD multi-slope models/DDexpMulti_geo_ADiag.R')
source('~/Dropbox/Temperate:Tropical/model_development/DD multi-slope models/DDlinMulti_geo_ADiag.R')
#source('~/Dropbox/Temperate:Tropical/analyses/Apr2018_wholeclade_approach/archived_troubleshooting scripts/16May_troubleshooting_for_Julien/createDDM_BETA.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/DDexp_geo_ADiag.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/DDlin_geo_ADiag.R', chdir = TRUE)

load('~/Dropbox/Manuscripts/PLOS Biology Tanager Project/data/stochastic maps for MCC tree/diet.simmaps.RData')

require(RPANDA)
library(corpcor)

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

o1E<-fit_t_comp_subgroup(smap,M,trim.class="fruit",regime.map=NULL,model="DDexp")
#> o1E
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
#[1] 3.094448
#
#$convergence
#[1] 0

o3E<-fit_t_DD(phylo=smap,data=M,model="exponential",subgroup="fruit",subgroup.map=smap,method="BB")
#> o3E
#$LogLik
#[1] -19.63071
#
#$AIC
#[1] 45.26141
#
#$AICc
#[1] 45.68999
#
#$rates
#            fruit       omn       inv
#beta  -0.01411994 0.0000000 0.0000000
#sigma  0.01105980 0.0110598 0.0110598
#
#$anc
#[1] 3.10385
#
#$convergence
#[1] 0
#
#$hess.values
#[1] 0
#
#$error
#NULL
#
#$param
#[1] -0.01411994 -4.50443869
#
#$phyloTrans
#NULL

o1L<-fit_t_comp_subgroup(smap,M,trim.class="fruit",regime.map=NULL,model="DDlin")
#> o1L
#$model
#[1] "DDlin"
#
#$LH
#[1] -19.77986
#
#$aic
#[1] 45.55971
#
#$aicc
#[1] 45.79272
#
#$free.parameters
#[1] 3
#
#$sig2
#[1] 0.009619716
#
#$r
#[1] -7.867985e-05
#
#$z0
#[1] 3.089411
#
#$convergence
#[1] 0

o3L<-fit_t_DD(phylo=smap,data=M,model="linear",subgroup="fruit",subgroup.map=smap,method="BB")
#> o3L
#$LogLik
#[1] -19.77547
#
#$AIC
#[1] 45.55095
#
#$AICc
#[1] 45.97952
#
#$rates
#              fruit         omn         inv
#beta  -7.851667e-05 0.000000000 0.000000000
#sigma  9.607082e-03 0.009607082 0.009607082
#
#$anc
#[1] 3.099127
#
#$convergence
#[1] 0
#
#$hess.values
#[1] 0
#
#$error
#NULL
#
#$param
#[1] -7.851667e-05  9.607082e-03
#
#$phyloTrans
#NULL


######################################
######################################
######################################
######################################

#two slope with subgroup pruning (see tests_20180929.R)

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
source('~/Dropbox/Scripts/R scripts/resortGeoObject.R', chdir = TRUE)
source('~/ddexp/tests/20180929/CreateGeoObject_simmap.subgroup.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/CreateSMatrix.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/ReconcileGeoObjectSMatrix.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/resortSMatrix.R', chdir = TRUE)
source('~/Dropbox/Temperate:Tropical/model_development/DD multi-slope models/DDexpMulti_geo_ADiag.R')
source('~/Dropbox/Temperate:Tropical/model_development/DD multi-slope models/DDlinMulti_geo_ADiag.R')

load('~/Dropbox/Manuscripts/PLOS Biology Tanager Project/data/stochastic maps for MCC tree/diet.simmaps.RData')

require(RPANDA)
library(corpcor)

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

o1<-fit_t_comp_subgroup(smap,M,trim.class="fruit",regime.map=smap2,model="DDexp")

##currently, can't get optim to search anything other than starting values, regardless of what they are

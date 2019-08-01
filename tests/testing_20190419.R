load("~/ddexp/data/motmot.data.Rd")
source('~/Dropbox/Scripts/R scripts/stratified_BGB_to_tables.R', chdir = TRUE)

source('~/ddexp/R/fit_t_EB.R', chdir = TRUE)
source('~/ddexp/R/common_functions.R', chdir = TRUE)
source('~/ddexp/R/generalized_functions.R', chdir = TRUE)
source('~/ddexp/R/fit_t_general.R', chdir = TRUE)
source('~/ddexp/R/fit_t_DD.R', chdir = TRUE)
source('~/ddexp/R/CreateClassbyClassObject_mvMORPH.r', chdir = TRUE)
source('~/ddexp/R/make.simmap.BGB.R', chdir = TRUE)
source('~/ddexp/R/CreateGeobyClassObject_mvMORPH.R', chdir = TRUE)
source('~/ddexp/R/CreateBioGeoB_Object_subclade.R', chdir = TRUE)
source('~/ddexp/R/stratified_BGB_to_tables.R')

smap<-stratified_BGB_to_tables(motmot.data[[1]],motmot.data[[2]],1)
geo.object<-CreateBioGeoB_Object_subclade(anc.phylo=motmot.data[[1]],subclade.phylo=motmot.data[[1]],ana.events=smap$ana.int,clado.events=smap$clado.int)


phylo<-motmot.data[[1]]
data<-motmot.data[[3]]
data<-data[which(!is.na(data))]
data<-data[which(names(data)%in%phylo$tip.label)]

phylo.trimmed<-keep.tip(phylo,names(data))

e1<-abs(rnorm(length(data),sd=0.05))
names(e1)<-phylo.trimmed$tip.label


source('~/Dropbox/Scripts/R scripts/DDexp_geo_ADiag_ME.R')
source('~/Dropbox/Scripts/R scripts/PhenotypicModel_PLUSME.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/PhenotypicADiag.R', chdir = TRUE)

ddexp_geo_me<-createModel_DDexp_geo_ME(tree=phylo.trimmed,geo.object=geo.object)

params<-c(0,log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))

fitTipData(ddexp_geo_me,data=data,error=e1,params0=params,GLSstyle=TRUE)

#*** Fit of tip trait data ***
#Finding the maximum likelihood estimator of the parameters, before returning the likelihood and the inferred parameters...
#Computation time : 0.716094 secs 
#$value
#[1] 6.469754
#
#$inferredParams
#         m0   logsigma0           r lognuisance 
#   4.354151   -2.263221   -0.247283   -2.145162 
#
#$convergence
#[1] 0

require(phytools)
require(geiger)
require(data.table)
geo.map<-make.simmap.BGB(anc.phylo=motmot.data[[1]],subclade.phylo=motmot.data[[1]],ana.events=smap$ana.int,clado.events=smap$clado.int,rnd=6)
fit_t_DD(phylo=phylo.trimmed,data=data,model="exponential",geo.map=geo.map[[1]],error=e1)
#
#$LogLik
#[1] -6.549749
#
#$AIC
#[1] 21.0995
#
#$AICc
#[1] 29.0995
#
#$rates
#                H          GH           G
#beta  0.024944964 0.024944964 0.024944964
#sigma 0.004669643 0.004669643 0.004669643
#
#$anc
#[1] 4.330262
#
#$convergence
#[1] 0
#
#$hess.values
#[1] 0
#
#$error
#[1] 0.07446634
#
#$param
#[1]  0.02494496 -5.36667273 -0.27288522
#
#$phyloTrans
#NULL

##Results are fairly similar (0.079995 lnL diff), but not identical, and slope is the opposite

##Still, likelihood is similar, so fit_t_DD isn't reaching ML that RPANDA is, but it's not far off:
> getDataLikelihood(ddexp_geo_me,data=data,error=e1,params=c(4.331188,log(sqrt(0.004619337)), 0.02628809, log(0.07294118)))
[1] 6.541768



##### TANAGER checks ####

load('~/ddexp/data/tanager.data.Rd')

require(RPANDA)
phylo<-tanager.data[[1]]
data<-tanager.data[[3]]
data<-data[which(!is.na(data))]
data<-data[which(names(data)%in%phylo$tip.label)]

phylo.trimmed<-keep.tip(phylo,names(data))

e1<-abs(rnorm(length(data),sd=0.05))
names(e1)<-phylo.trimmed$tip.label

geography.object<-CreateGeoObject_BioGeoBEARS(full.phylo=phylo,trimmed.phylo=phylo.trimmed,ana.events=tanager.data[[2]][[1]][[1]],clado.events=tanager.data[[2]][[2]][[1]])

source('~/Dropbox/Scripts/R scripts/DDexp_geo_ADiag_ME.R')

ddexp_geo_me<-createModel_DDexp_geo_ME(tree=phylo.trimmed,geo.object=geography.object)

params<-c(0,log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))

fitTipData(ddexp_geo_me,data=data,error=e1,params0=params,GLSstyle=TRUE)
#*** Fit of tip trait data ***
#Finding the maximum likelihood estimator of the parameters, before returning the likelihood and the inferred parameters...
#Computation time : 2.189974 hours 
#$value
#[1] 109.6904
#
#$inferredParams
#         m0   logsigma0           r lognuisance 
# 3.03826367 -2.38642860 -0.01115995 -3.55587031 
#
#$convergence
#[1] 0


source('~/ddexp/R/make.simmap.BGB.R')
require(phytools)
require(geiger)
require(data.table)
geo.map<-make.simmap.BGB(anc.phylo=phylo,subclade.phylo=phylo.trimmed,ana.events=tanager.data[[2]][[1]][[1]],clado.events=tanager.data[[2]][[2]][[1]],rnd=6)

fit_t_DD(phylo=phylo.trimmed,data=data,model="exponential",geo.map=geo.map[[1]],error=e1)

#$LogLik
#[1] -109.676
#
#$AIC
#[1] 227.3521
#
#$AICc
#[1] 227.4782
#
#$rates
#              EIFC          IFC            I           ID           IG
#beta  -0.005924849 -0.005924849 -0.005924849 -0.005924849 -0.005924849
#sigma  0.006591190  0.006591190  0.006591190  0.006591190  0.006591190
#
#$anc
#[1] 3.043738
#
#$convergence
#[1] 0
#
#$hess.values
#[1] 0
#
#$error
#[1] 0.02507703
#
#$param
#[1] -0.005924849 -5.022021433  0.158357283
#
#$phyloTrans
#NULL

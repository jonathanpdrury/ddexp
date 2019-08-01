set.seed(1224)
t1<-rcoal(50)
d1<-rnorm(50)
names(d1)<-t1$tip.label
e1<-rnorm(50,sd=0.05)
names(e1)<-t1$tip.label


source('~/Dropbox/Scripts/R scripts/PhenotypicModel_PLUSME.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/PhenotypicADiag.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/DDexp_nogeo_ADiag_ME.R', chdir = TRUE)

ddexp_me<-createModel_DDexp_ME(t1)
fitTipData(ddexp_me,data=d1,error=e1,GLSstyle=TRUE)

#*** Fit of tip trait data ***
#Finding the maximum likelihood estimator of the parameters, before returning the likelihood and the inferred parameters...
#Computation time : 4.571577 secs 
#$value
#[1] 71.9436
#
#$inferredParams
#         m0   logsigma0           r lognuisance 
# 0.09313754 -4.75617415 -3.98344760  0.03764250 
#
#$convergence
#[1] 0


source('~/ddexp/R/fit_t_DD.R')
source('~/ddexp/R/fit_t_general.R')
source('~/ddexp/R/generalized_functions.R')

fit_t_DD(phylo=t1,data=d1,model="exponential",error=e1)
#$LogLik
#[1] -71.9436
#
#$AIC
#[1] 151.8872
#
#$AICc
#[1] 152.7761
#
#$rates
#               [,1]
#beta  -2.710479e+00
#sigma  7.849609e-06
#
#$anc
#[1] 0.09313745
#
#$convergence
#[1] 0
#
#$hess.values
#[1] 0
#
#$error
#[1] 1.038019
#
#$param
#[1]  -2.710479 -11.755047   1.018832
#
#$phyloTrans
#NULL



### Now check DD with geo scripts
### NOTE, neither approach is working right now, should trim to something smaller for testing

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
#Computation time : 0.6200631 secs 
#$value
#[1] -0.8226077
#
#$inferredParams
#         m0   logsigma0           r lognuisance 
# -0.4966965  -1.2870673  -0.7589808 -24.3747843 
#
#$convergence
#[1] 0

source('~/ddexp/R/make.simmap.BGB.R')
require(phytools)
require(geiger)
require(data.table)
geo.map<-make.simmap.BGB(anc.phylo=phylo,subclade.phylo=phylo.trimmed,ana.events=tanager.data[[2]][[1]][[1]],clado.events=tanager.data[[2]][[2]][[1]],rnd=6)
fit_t_DD(phylo=phylo.trimmed,data=pPC1,model="exponential",geo.map=geo.map,error=error)
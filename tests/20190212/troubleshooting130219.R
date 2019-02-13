##these tests check same fits with wrapper script, which is now working by passing class.df and class.objects directly to fun[[funInd]](x,df,times)


######################################
######################################
######################################
######################################

##testing two slope version

######################################
######################################
######################################
######################################

require(phytools)
require(RPANDA)

set.seed(4321)

fiftytipyuletrees<-pbtree(n=50,nsim=1)
regions<-c("trop","trop","trop","temp","temp")
phylo<-fiftytipyuletrees
samp.reg<-sample(regions)
trait<-c(rep(samp.reg[1],10),rep(samp.reg[2],10),rep(samp.reg[3],10),rep(samp.reg[4],10),rep(samp.reg[5],10))
names(trait)<-phylo$tip.label
smap<-make.simmap(phylo,trait,model="ARD")
source('~/Dropbox/Temperate:Tropical/analyses/model_development/DD multi-slope models/sim_DDmulti.R')



source('~/Dropbox/Temperate:Tropical/analyses/model_development/DD multi-slope models/DDexpMulti_nogeo_ADiag.R')
source('~/Dropbox/Temperate:Tropical/analyses/model_development/DD multi-slope models/DDlinMulti_nogeo_ADiag.R')

source('~/Dropbox/Scripts/R scripts/resortSMatrix.R')
source('~/Dropbox/Scripts/R scripts/CreateSMatrix.R')

r.matrix= CreateSMatrix(CreateClassObject(smap),S.cats=c("trop","temp"))
trait<-try(sim_multiDD(smap,pars=matrix(c(0.05,-0.05,-0.1,0),nrow=1),plot=F,S.matrix=r.matrix,rnd=5,model="exponential"))

pc<-proc.time()
rmats<-resortSMatrix(phylo,r.matrix)
ddm.ob<-createModel_DDexp_multi(phylo,rmats)
params0<-c(0,log(sqrt(var(trait[1,])/max(nodeHeights(phylo)))),-1e-4,-1e-4)	
o1a<-fitTipData(ddm.ob,trait[1,],params0,GLSstyle=TRUE)


ddl.ob<-createModel_DDlin_multi(phylo,rmats)
o1al<-fitTipData(ddl.ob,trait[1,],params0,GLSstyle=TRUE)
print(proc.time()-pc)
#> print(proc.time()-pc)
#   user  system elapsed 
#  8.901   0.148   9.002 

#$value
#[1] -42.75546
#
#$inferredParams
#          m0    logsigma0           r1           r2 
#-0.009016524 -1.724855606 -0.046974344 -0.071106515 
#
#$convergence
#[1] 0


#$value
#[1] -42.63477
#
#$inferredParams
#           m0     logsigma0            r1            r2 
#-0.0066771301 -1.8346771805 -0.0006239160 -0.0008672482 
#
#$convergence
#[1] 0

source('~/ddexp/R/fit_t_DD_BETA.R', chdir = TRUE)
source('~/ddexp/R/fit_t_general_options_old_JPD_130219.R', chdir = TRUE)
source('~/ddexp/R/generalized_functions.R')
source('~/Dropbox/Scripts/R scripts/CreateClassObject.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/trimSimmap.R', chdir = TRUE)
source('~/ddexp/R/CreateClassbyClassObject_mvMORPH.R')


pc<-proc.time()
fit_t_DD(phylo=phylo,data=trait[1,],model="both",regime.map=smap,method="Nelder-Mead")
print(proc.time()-pc)
#> print(proc.time()-pc)
#   user  system elapsed 
#  4.307   0.060   4.316 

#$exponential.fit
#$exponential.fit$LogLik
#[1] 42.76091
#
#$exponential.fit$AIC
#[1] -77.52181
#
#$exponential.fit$AICc
#[1] -76.63293
#
#$exponential.fit$rates
#             temp        trop
#beta  -0.07119071 -0.04711833
#sigma  0.03182148  0.03182148
#
#$exponential.fit$anc
#[1] -0.009018124
#
#$exponential.fit$convergence
#[1] 0
#
#$exponential.fit$hess.values
#[1] 0
#
#$exponential.fit$error
#NULL
#
#$exponential.fit$param
#[1] -0.07119071 -0.04711833 -3.44761375
#
#$exponential.fit$phyloTrans
#NULL
#
#
#$linear.fit
#$linear.fit$LogLik
#[1] 42.64289
#
#$linear.fit$AIC
#[1] -77.28578
#
#$linear.fit$AICc
#[1] -76.39689
#
#$linear.fit$rates
#               temp          trop
#beta  -0.0008692097 -0.0006261159
#sigma  0.0255324914  0.0255324914
#
#$linear.fit$anc
#[1] -0.00665851
#
#$linear.fit$convergence
#[1] 0
#
#$linear.fit$hess.values
#[1] 0
#
#$linear.fit$error
#NULL
#
#$linear.fit$param
#[1] -0.0008692097 -0.0006261159 -3.6678034642
#
#$linear.fit$phyloTrans
#NULL


###two-slope version is working fine for both cases

######################################
######################################
######################################
######################################

##one slope version
##not working


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
set.seed(1234)
data<-rnorm(length(BGB.examples$Canidae.phylo$tip.label))
names(data)<-BGB.examples$Canidae.phylo$tip.label

pc<-proc.time()
DDexp.geo.fit<-fit_t_comp(BGB.examples$Canidae.phylo, data, model="DDexp", geography.object=Canidae.geography.object)
DDlin.geo.fit<-fit_t_comp(BGB.examples$Canidae.phylo, data,model="DDlin", geography.object=Canidae.geography.object)
print(proc.time()-pc)
#   user  system elapsed 
#  3.573   0.051   3.600 

#> DDexp.geo.fit 
#$LH
#[1] -61.96697
#
#$aic
#[1] 129.9339
#
#$aicc
#[1] 130.7339
#
#$free.parameters
#[1] 3
#
#$sig2
#[1] 0.1123035
#
#$r
#[1] 0.2169455
#
#$z0
#[1] -0.2401745
#
#$convergence
#[1] 0
#
#> DDlin.geo.fit
#$LH
#[1] -61.26467
#
#$aic
#[1] 128.5293
#
#$aicc
#[1] 129.3293
#
#$free.parameters
#[1] 3
#
#$sig2
#[1] 6.624159e-06
#
#$b
#[1] 0.1100047
#
#$z0
#[1] -0.2475576
#
#$convergence
#[1] 0


###checking fit_t_DD
source('~/ddexp/R/fit_t_DD_BETA.R', chdir = TRUE)
source('~/ddexp/R/make.simmap.BGB.R', chdir = TRUE)
source('~/ddexp/R/fit_t_general_options_old_JPD_130219.R', chdir = TRUE)
source('~/ddexp/R/generalized_functions.R')


require(geiger)
require(data.table)
require(phytools)
geo.map=make.simmap.BGB(BGB.examples$Canidae.phylo,BGB.examples$Canidae.phylo,ana.events=BGB.examples$Canidae.ana.events,clado.events=BGB.examples$Canidae.clado.events,rnd=6)$geo.simmap

pc<-proc.time()
fit_t_DD(phylo,data=data,model="both",geo.map=geo.map,method=c("Nelder-Mead"))
print(proc.time()-pc)
#   user  system elapsed 
#  1.929   0.031   1.952 
#
#$exponential.fit
#$exponential.fit$LogLik
#[1] -61.96749
#
#$exponential.fit$AIC
#[1] 129.935
#
#$exponential.fit$AICc
#[1] 130.735
#
#$exponential.fit$rates
#              F        DF         D        BD       BCD      BCDH       BDH     BCDGH      CDGH         C        FI
#beta  0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598
#sigma 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631
#              I       FEI        BI       ABI      ABCI     ABCDI    ABCDFI        CD       ACD      ABCD        CH
#beta  0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598
#sigma 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631
#              B        AB       BFG        BJ        FG         G        BF       BIJ        GH       FGH        FH
#beta  0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598
#sigma 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631
#            CFH      CDFH        AD       CDH       FIJ      AFIJ
#beta  0.2170598 0.2170598 0.2170598 0.2170598 0.2170598 0.2170598
#sigma 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631 0.1121631
#
#$exponential.fit$anc
#[1] -0.2400122
#
#$exponential.fit$convergence
#[1] 0
#
#$exponential.fit$hess.values
#[1] 0
#
#$exponential.fit$error
#NULL
#
#$exponential.fit$param
#[1]  0.2170598 -2.1878009
#
#$exponential.fit$phyloTrans
#NULL
#
#
#$linear.fit
#$linear.fit$LogLik
#[1] -61.26314
#
#$linear.fit$AIC
#[1] 128.5263
#
#$linear.fit$AICc
#[1] 129.3263
#
#$linear.fit$rates
#                 F           DF            D           BD          BCD         BCDH          BDH        BCDGH
#beta  1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01
#sigma 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08
#              CDGH            C           FI            I          FEI           BI          ABI         ABCI
#beta  1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01
#sigma 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08
#             ABCDI       ABCDFI           CD          ACD         ABCD           CH            B           AB
#beta  1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01
#sigma 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08
#               BFG           BJ           FG            G           BF          BIJ           GH          FGH
#beta  1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01
#sigma 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08
#                FH          CFH         CDFH           AD          CDH          FIJ         AFIJ
#beta  1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01 1.101020e-01
#sigma 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08 7.481076e-08
#
#$linear.fit$anc
#[1] -0.2471769
#
#$linear.fit$convergence
#[1] 0
#
#$linear.fit$hess.values
#[1] 0
#
#$linear.fit$error
#NULL
#
#$linear.fit$param
#[1]   0.110102 -16.408304
#
#$linear.fit$phyloTrans
#NULL


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
source('~/ddexp/R/stratified_BGB_to_tables.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/resortGeoObject.R', chdir = TRUE)
source('~/ddexp/tests/20180929/CreateGeoObject_simmap.subgroup.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/CreateSMatrix.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/ReconcileGeoObjectSMatrix.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/resortSMatrix.R', chdir = TRUE)
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

pc<-proc.time()
o1E<-fit_t_comp_subgroup(smap,M,trim.class="fruit",regime.map=NULL,model="DDexp")
o1L<-fit_t_comp_subgroup(smap,M,trim.class="fruit",regime.map=NULL,model="DDlin")
print(proc.time()-pc)
#   user  system elapsed 
#183.705  12.971 197.108 

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

source('~/ddexp/R/fit_t_DD_BETA.R', chdir = TRUE)
source('~/ddexp/R/fit_t_general_options_old_JPD_130219.R', chdir = TRUE)
source('~/ddexp/R/generalized_functions.R')

pc<-proc.time()
fit_t_DD(phylo=smap,data=M,model="both",subgroup="fruit",subgroup.map=smap,method="Nelder-Mead")
print(proc.time()-pc)
#   user  system elapsed 
#  3.260   0.077   3.339 

#$exponential.fit
#$exponential.fit$LogLik
#[1] -19.63071
#
#$exponential.fit$AIC
#[1] 45.26141
#
#$exponential.fit$AICc
#[1] 45.68999
#
#$exponential.fit$rates
#            fruit        omn        inv
#beta  -0.01412492 0.00000000 0.00000000
#sigma  0.01106355 0.01106355 0.01106355
#
#$exponential.fit$anc
#[1] 3.103856
#
#$exponential.fit$convergence
#[1] 0
#
#$exponential.fit$hess.values
#[1] 0
#
#$exponential.fit$error
#NULL
#
#$exponential.fit$param
#[1] -0.01412492 -4.50409901
#
#$exponential.fit$phyloTrans
#NULL
#
#
#$linear.fit
#$linear.fit$LogLik
#[1] -19.77547
#
#$linear.fit$AIC
#[1] 45.55094
#
#$linear.fit$AICc
#[1] 45.97952
#
#$linear.fit$rates
#              fruit         omn         inv
#beta  -7.860514e-05 0.000000000 0.000000000
#sigma  9.612849e-03 0.009612849 0.009612849
#
#$linear.fit$anc
#[1] 3.099135
#
#$linear.fit$convergence
#[1] 0
#
#$linear.fit$hess.values
#[1] 0
#
#$linear.fit$error
#NULL
#
#$linear.fit$param
#[1] -7.860514e-05 -4.644655e+00
#
#$linear.fit$phyloTrans
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
source('~/ddexp/R/make.simmap.BGB.R', chdir = TRUE)
source('~/ddexp/R/stratified_BGB_to_tables.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/resortGeoObject.R', chdir = TRUE)
source('~/ddexp/tests/20180929/CreateGeoObject_simmap.subgroup.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/CreateSMatrix.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/ReconcileGeoObjectSMatrix.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/resortSMatrix.R', chdir = TRUE)
source('~/Dropbox/Temperate:Tropical/analyses/model_development/DD multi-slope models/DDexpMulti_geo_ADiag.R')
source('~/Dropbox/Temperate:Tropical/analyses/model_development/DD multi-slope models/DDlinMulti_geo_ADiag.R')

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

pc<-proc.time()
o1e<-fit_t_comp_subgroup(smap,M,trim.class="fruit",regime.map=smap2,model="DDexp")
o1l<-fit_t_comp_subgroup(smap,M,trim.class="fruit",regime.map=smap2,model="DDlin")
print(proc.time()-pc)
#   user  system elapsed 
#366.061  16.301 383.818 

#$model
#[1] "DDexp"
#
#$LH
#[1] -19.55769
#
#$aic
#[1] 45.11538
#
#$aicc
#[1] 45.34839
#
#$free.parameters
#[1] 3
#
#$sig2
#[1] 0.01114394
#
#$r1
#[1] -0.02652755
#
#$r2
#[1] -0.0266374
#
#$z0
#[1] 3.095275
#
#$convergence
#[1] 0


#$model
#[1] "DDlin"
#
#$LH
#[1] -20.46297
#
#$aic
#[1] 46.92595
#
#$aicc
#[1] 47.15896
#
#$free.parameters
#[1] 3
#
#$sig2
#[1] 0.005889312
#
#$r1
#[1] 2.861665e-05
#
#$r2
#[1] 2.31407e-07
#
#$z0
#[1] 3.089581
#
#$convergence
#[1] 0

source('~/ddexp/R/fit_t_DD_BETA.R', chdir = TRUE)
source('~/ddexp/R/fit_t_general_options_old_JPD_130219.R', chdir = TRUE)
source('~/ddexp/R/generalized_functions.R')

pc<-proc.time()
fit_t_DD(phylo=smap,data=M,model="both",subgroup="fruit",subgroup.map=smap,regime.map=smap2,method="Nelder-Mead")
print(proc.time()-pc)
#   user  system elapsed 
#  7.459   0.128   7.586 


#$exponential.fit
#$exponential.fit$LogLik
#[1] -19.53093
#
#$exponential.fit$AIC
#[1] 47.06187
#
#$exponential.fit$AICc
#[1] 47.78914
#
#$exponential.fit$rates
#        temperate          Z    tropical
#beta  -0.02649129 0.00000000 -0.02673028
#sigma  0.01113273 0.01113273  0.01113273
#
#$exponential.fit$anc
#[1] 3.096731
#
#$exponential.fit$convergence
#[1] 0
#
#$exponential.fit$hess.values
#[1] 0
#
#$exponential.fit$error
#NULL
#
#$exponential.fit$param
#[1] -0.02649129 -0.02673028 -4.49786621
#
#$exponential.fit$phyloTrans
#NULL
#
#
#$linear.fit
#$linear.fit$LogLik
#[1] -20.46343
#
#$linear.fit$AIC
#[1] 48.92686
#
#$linear.fit$AICc
#[1] 49.65413
#
#$linear.fit$rates
#         temperate           Z     tropical
#beta  2.834454e-05 0.000000000 2.297838e-07
#sigma 5.891906e-03 0.005891906 5.891906e-03
#
#$linear.fit$anc
#[1] 3.089635
#
#$linear.fit$convergence
#[1] 0
#
#$linear.fit$hess.values
#[1] 0
#
#$linear.fit$error
#NULL
#
#$linear.fit$param
#[1]  2.834454e-05  2.297838e-07 -5.134176e+00
#
#$linear.fit$phyloTrans
#NULL

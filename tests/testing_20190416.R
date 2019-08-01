
source('~/Dropbox/Scripts/R scripts/PhenotypicModel_JPDv2.R', chdir = TRUE)
#source('~/Dropbox/Scripts/R scripts/PhenotypicModel_PLUSME.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/PhenotypicADiag.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/DDexp_nogeo_ADiag.R', chdir = TRUE)
#source('~/Dropbox/Scripts/R scripts/DDexp_nogeo_ADiag_ME.R', chdir = TRUE)

#TESTING VERSION WHERE ERROR ISNT ADDED
#            V<- V.int #+ diag(error^2) + diag(rep(exp(params[length(params)]),n))

require(ape)
require(phytools)

set.seed(1234)
t1<-rcoal(3)

set.seed(1234)
d1<-rnorm(3)
names(d1)<-t1$tip.label

set.seed(1234)
e1<-rnorm(3,sd=0.04)
e1<-abs(e1)
names(e1)<-t1$tip.label

m1<-createModel_DDexp(t1)
#m1<-createModel_DDexp_ME(t1)

params=c(0,0,0)
#params=c(0,0,0,0)

getTipDistribution(m1,params)
#$mean
#   [,1]
#t3    0
#t1    0
#t2    0
#
#$Sigma
#          t3       t1        t2
#t3 1.0806784 0.000000 0.2467589
#t1 0.0000000 1.080678 0.0000000
#t2 0.2467589 0.000000 1.0806784

params2=c(0,0,-1)
#params2=c(0,0,-1,0)
getTipDistribution(m1,params2)
#$mean
#   [,1]
#t3    0
#t1    0
#t2    0
#
#$Sigma
#           t3         t1         t2
#t3 0.07491359 0.00000000 0.03339518
#t1 0.00000000 0.07491359 0.00000000
#t2 0.03339518 0.00000000 0.07491359

#params3=c(0,0,-1,1)
#getTipDistribution(m1,params3)

params4=c(0,log(0.05),-0.01)

#params4=c(0,log(0.05),-0.01,0)

getTipDistribution(m1,params4)
#$mean
#   [,1]
#t3    0
#t1    0
#t2    0
#
#$Sigma
#             t3          t1           t2
#t3 0.0026278655 0.000000000 0.0006046818
#t1 0.0000000000 0.002627866 0.0000000000
#t2 0.0006046818 0.000000000 0.0026278655

params5=c(0.05154108, -23.39056447,  15.61937201)

getTipDistribution(m1,params5)

getDataLikelihood(m1,d1,e1,params)
#[1] 4.213949

getDataLikelihood(m1,d1,e1,params2)
#[1] 21.87242

getDataLikelihood(m1,d1,e1,params3)
#[1] 21.87242

getDataLikelihood(m1,d1,e1,params4)
#[1] 556.7222

getDataLikelihood(m1,d1,e1,params5)
#[1] 4.099885


###Likelihoods and tip distributions are exactly the same, but the fitTipDistribution output isn't, though they do return exact likelihood (on ridge?).

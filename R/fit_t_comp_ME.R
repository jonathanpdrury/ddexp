source('PhenotypicModel_PLUSME.R')
source('PhenotypicADiag.R')
source('DDlinMulti_geo_ADiag_ME.R')
source('DDlinMulti_nogeo_ADiag_ME.R')
source('DDexpMulti_nogeo_ADiag_ME.R')
source('DDexpMulti_geo_ADiag_ME.R')
source('MC_twoS_PM_geo_ME.R')
source('MC_twoS_PM_ME.R')
source('DDexp_geo_ADiag_ME.R')
source('DDlin_geo_ADiag_ME.R')
source('MC_geo_PM_ME.R')
source('DDexp_nogeo_ADiag_ME.R')
source('MC_nogeo_ADiag_ME.R')
source('DDlin_nogeo_ADiag_ME.R')
source('resortGeoObject.R')
source('resortSMatrix.R')

#source('~/Dropbox/Scripts/R scripts/PhenotypicModel_PLUSME.R')
#source('~/Dropbox/Scripts/R scripts/PhenotypicADiag.R')
#source('~/Dropbox/Scripts/R scripts/DDlinMulti_geo_ADiag_ME.R')
#source('~/Dropbox/Scripts/R scripts/DDlinMulti_nogeo_ADiag_ME.R')
#source('~/Dropbox/Scripts/R scripts/DDexpMulti_nogeo_ADiag_ME.R')
#source('~/Dropbox/Scripts/R scripts/DDexpMulti_geo_ADiag_ME.R')
#source('~/Dropbox/Scripts/R scripts/MC_twoS_PM_geo_ME.R')
#source('~/Dropbox/Scripts/R scripts/MC_twoS_PM_ME.R')
#source('~/Dropbox/Scripts/R scripts/DDexp_geo_ADiag_ME.R')
#source('~/Dropbox/Scripts/R scripts/DDlin_geo_ADiag_ME.R')
#source('~/Dropbox/Scripts/R scripts/MC_geo_PM_ME.R')
#source('~/Dropbox/Scripts/R scripts/DDexp_nogeo_ADiag_ME.R')
#source('~/Dropbox/Scripts/R scripts/MC_nogeo_ADiag_ME.R')
#source('~/Dropbox/Scripts/R scripts/DDlin_nogeo_ADiag_ME.R')
#source('~/Dropbox/Scripts/R scripts/resortGeoObject.R')
#source('~/Dropbox/Scripts/R scripts/resortSMatrix.R')

fit_t_comp_ME<-function(phylo,data,error, model=c("MC","DDexp","DDlin"),pars=NULL,geography.object=NULL, regime.map=NULL){

#check to make sure data are univariate, with names matching phylo object
if(length(data)!=length(phylo$tip.label)){stop("length of data does not match length of tree")}
if(is.null(names(data))){stop("data missing taxa names")}
if(!is.null(dim(data))){stop("data needs to be a single trait")}
is_tip <- phylo$edge[,2] <= length(phylo$tip.label)
if(sum(diff(phylo$edge[is_tip, 2])<0)>0){ stop('fit_t_comp cannot be used with ladderized phylogenies')}


if(is.null(geography.object) & is.null(regime.map)){ #single-slope version for sympatric clades


	if(model=="MC"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),-0.1,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		mc.ob<-createModel_MC_ME(phylo)
		opt<-fitTipData(mc.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		S<-opt$inferredParams[3]
		nuisance<-exp(opt$inferredParams[4])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*4 - 2*(-opt$value)), aicc = (2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), S = as.numeric(S), nuisance=as.numeric(nuisance),z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddexp.ob<-createModel_DDexp_ME(phylo)
		opt<-fitTipData(ddexp.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		r<-opt$inferredParams[3]
		nuisance<-exp(opt$inferredParams[4])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), r = as.numeric(r), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddlin.ob<-createModel_DDlin_ME(phylo)
		opt<-fitTipData(ddlin.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		b<-opt$inferredParams[3]
		nuisance<-exp(opt$inferredParams[4])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), b = as.numeric(b), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
}

if(!is.null(geography.object) & is.null(regime.map)){ #single-slope version with biogeography


	if(length(geography.object$geography.object)<phylo$Nnode){stop("geography object cannot have fewer components than internode intervals in phylo")}
	sgeo<-resortGeoObject(phylo,geography.object) #resorts geo.object to match tip label order in Marc code
	if(model=="MC"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),-0.1,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		mc.ob<-createModel_MC_geo_ME(phylo,sgeo)
		opt<-fitTipData(mc.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		S<-opt$inferredParams[3]
		nuisance<-exp(opt$inferredParams[4])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc =(2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), S = as.numeric(S), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddexp.ob<-createModel_DDexp_geo_ME(phylo,sgeo)
		opt<-fitTipData(ddexp.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		r<-opt$inferredParams[3]
		nuisance<-exp(opt$inferredParams[4])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), r = as.numeric(r), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddlin.ob<-createModel_DDlin_geo_ME(phylo,sgeo)
		opt<-fitTipData(ddlin.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		b<-opt$inferredParams[3]
		nuisance<-exp(opt$inferredParams[4])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*4 - 2*(-opt$value))+((2*4*(4+1))/(length(phylo$tip.label)-4-1)), free.parameters = 4, sig2 = as.numeric(sig2), b = as.numeric(b), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}

}

if(is.null(geography.object) & !is.null(regime.map)){ #multi-slope version for sympatric clades

	smat<-resortSMatrix(phylo,regime.map)
	if(model=="MC"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),-0.1,-0.1,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		mc.ob<-createModel_MC_twoS_ME(phylo,smat)
		opt<-fitTipData(mc.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		S1<-opt$inferredParams[3]
		S2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)), free.parameters = 5, sig2 = as.numeric(sig2), S1 = as.numeric(S1), S2 = as.numeric(S2), nuisance=as.numeric(nuisance),z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddexp.ob<-createModel_DDexp_multi_ME(phylo,smat)
		opt<-fitTipData(ddexp.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		r1<-opt$inferredParams[3]
		r2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)), free.parameters = 5, sig2 = as.numeric(sig2), r1 = as.numeric(r1), r2 = as.numeric(r2), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddlin.ob<-createModel_DDlin_multi_ME(phylo,smat)
		opt<-fitTipData(ddlin.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		b1<-opt$inferredParams[3]
		b2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)), free.parameters = 5, sig2 = as.numeric(sig2), b1 = as.numeric(b1),b2 = as.numeric(b2), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
}

if(!is.null(geography.object) & !is.null(regime.map)){ #multi-slope version with biogeography

if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),-0.1,-0.1,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
params0<-c(0,pars)

	if(length(geography.object$geography.object)<phylo$Nnode){stop("geography object cannot have fewer components than internode intervals in phylo")}
	sgeo<-resortGeoObject(phylo,geography.object) #resorts geo.object to match tip label order in Marc code
	smat<-resortSMatrix(phylo,regime.map)
	if(model=="MC"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),-0.1,-0.1,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		mc.ob<-createModel_MC_twoS_geo_ME(phylo,sgeo,smat)
		opt<-fitTipData(mc.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		S1<-opt$inferredParams[3]
		S2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)), free.parameters = 5, sig2 = as.numeric(sig2), S1 = as.numeric(S1), S2 = as.numeric(S2), nuisance=as.numeric(nuisance),z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDexp"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddexp.ob<-createModel_DDexp_multi_geo_ME(phylo,sgeo,smat)
		opt<-fitTipData(ddexp.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		r1<-opt$inferredParams[3]
		r2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)), free.parameters = 5, sig2 = as.numeric(sig2), r1 = as.numeric(r1), r2 = as.numeric(r2), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}
	if(model=="DDlin"){
		if(is.null(pars)){pars<-c(log(0.95*sqrt(var(data)/max(nodeHeights(phylo)))),0,0,log(0.05*sqrt(var(data)/max(nodeHeights(phylo)))))}
		params0<-c(0,pars)
		ddlin.ob<-createModel_DDlin_multi_geo_ME(phylo,sgeo,smat)
		opt<-fitTipData(ddlin.ob,data,error,params0=params0,GLSstyle=TRUE)
		sig2<-(exp(opt$inferredParams[2]))^2
		b1<-opt$inferredParams[3]
		b2<-opt$inferredParams[4]
		nuisance<-exp(opt$inferredParams[5])
		z0<-opt$inferredParams[1]
		results<-list(LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*5 - 2*(-opt$value))+((2*5*(5+1))/(length(phylo$tip.label)-5-1)), free.parameters = 5, sig2 = as.numeric(sig2), b1 = as.numeric(b1),b2 = as.numeric(b2), nuisance=as.numeric(nuisance), z0 = as.numeric(z0), convergence = opt$convergence)
		return(results)
		}

}

}
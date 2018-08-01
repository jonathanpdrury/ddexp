### THIS VERSION FITS data to first stochastic map of both biogeography and of guild;
### final publication versions should parallelize/loop over combinations

load('tanagers_RES_ana_events_tables_SMv2.RData')
load('tanagers_RES_clado_events_tables_SMv2.RData')
load('diet.simmaps.RData')

require(deSolve)
require(geiger)
require(phytools)
library(corpcor)
require(methods)
require(tools)

source('CreateBioGeoB_Object_subclade_V2.R')
source('PhenotypicModel_JPDv2.R')
source('PhenotypicADiag.R')
source('DDexp_geo_ADiag.R')
source('resortGeoObject.R')
source('trimSimmap.R')
source('CreateClassObject.R')
source('CreateGeobyClassObject.R')
source('likelihood_subguild.R')
source('fit_t_comp_subguild.R')

source('fit_t_general.R')
source('generalized_functions.R')
source('CreateGeobyClassObject_mvMORPH.R')
source('make.simmap.BGB.R')

tan.tree<-read.tree("tanagertree_10regions.tree")
tandata<-read.csv("MASTER_tandata.csv")

guild.map<-diet.simmaps


fruit<-subset(tandata, diet=="fruit")#62 spp, #59 in tree
mass<-fruit[,5]
names(mass)<-fruit[,3]
nc<-name.check(tan.tree,mass)
fruit.tree<-drop.tip(tan.tree,nc$tree_not_data)
subdata<-fruit[which(fruit$treename%in%fruit.tree$tip.label),]


N=2
j=5:11
res.mat<-matrix(nrow=length(j)*N,ncol=20)
colnames(res.mat)<-c("subset","geo.map","trait","N","BM.log_lik","BM.sig2","BM.z0","DDexpgeo.iter","DDexpgeo_RPANDA.lnL","DDexpgeo_RPANDA.sig2","DDexpgeo_RPANDA.r","DDexpgeo_RPANDA.z0","DDexpgeo_RPANDA.conv","DDexpgeo_RPANDA.time","DDexpgeo_mvMORPH.lnL","DDexpgeo_mvMORPH.sig2","DDexpgeo_mvMORPH.r","DDexpgeo_mvMORPH.z0","DDexpgeo_mvMORPH.conv","DDexpgeo_mvMORPH.time")
counter=1

for(i in 1:2){
	guild.mapi<-drop.tip.simmap(guild.map[[i]],guild.map[[i]]$tip.label[which(!guild.map[[i]]$tip.label%in%tan.tree$tip.label)])
	for(j in 5:11){
		iter=i	
		M<-subdata[,j]
		subtree<-fruit.tree
		names(M)<-subdata$treename
		M<-subset(M,M!='NA')

		nc<-name.check(subtree,M)
		nc
		if(is.list(nc)){
		subtree<-drop.tip(subtree,nc$tree_not_data)
		}

		o2<-fitContinuous(subtree,M,model="BM")
		BM.log_lik<-o2$opt$lnL
		BM.sig2<-o2$opt$sigsq
		BM.z0<-o2$opt$z0
				
		iter=i
		guild.mapi<-drop.tip.simmap(guild.map[[i]],guild.map[[i]]$tip.label[which(!guild.map[[i]]$tip.label%in%tan.tree$tip.label)])
		pt<-proc.time()
		o5<-try(fit_t_comp_subguild(full.phylo=tan.tree,ana.events=tanagers_RES_ana_events_tables_SMv2[[i]],clado.events=tanagers_RES_clado_events_tables_SMv2[[i]],map=guild.mapi,data=M,trim.class="fruit",model="DDexp",pars=NULL,method="Nelder-Mead",bounds=NULL))
		while(class(o5)=="try-error"){
			iter=iter+1
			guild.mapi<-drop.tip.simmap(guild.map[[iter]],guild.map[[iter]]$tip.label[which(!guild.map[[iter]]$tip.label%in%tan.tree$tip.label)])
			o5<-try(fit_t_comp_subguild(full.phylo=tan.tree,ana.events=tanagers_RES_ana_events_tables_SMv2[[i]],clado.events=tanagers_RES_clado_events_tables_SMv2[[i]],map=guild.mapi,data=M,trim.class="fruit",model="DDexp",pars=NULL,method="Nelder-Mead",bounds=NULL))
			}
		DDexpgeo.iter<-iter	
		DDexpgeo_RPANDA.lnL<-o5$LH
		DDexpgeo_RPANDA.sig2<-o5$sig2
		DDexpgeo_RPANDA.r<-o5$r
		DDexpgeo_RPANDA.z0<-o5$z0
		DDexpgeo_RPANDA.conv<-o5$convergence
		
		DDexpgeo_RPANDA.time<-proc.time()-pt
		
		###now add in fit_t_* tools and record output as well as time
		out<-CreateGeobyClassObject_mvMORPH(phylo=tan.tree,map=guild.mapi,ana.events=tanagers_RES_ana_events_tables_SMv2[[i]],clado.events=tanagers_RES_clado_events_tables_SMv2[[i]],trim.class="fruit")
		geo.class.df<-return.class.df(out$geo.simmap,out$geo.class.object)
		geo.class.df[,which(colnames(out$geo.simmap$mapped.edge)=='Z')+1]=0
		##to check against BM, can use this approach: 		geo.class.df[,2:dim(geo.class.df)[2]][!is.na(geo.class.df[,2:dim(geo.class.df)[2]])]=0

		new_list_function<-create.function.list(out$geo.simmap,out$geo.class.object,geo.class.df)

		pt<-proc.time()
		o3<-fit_t_general_subgroup(out$geo.simmap, M, fun=new_list_function,diagnostic=T,echo=T)
##NOTE: this is where the problem was, tips shouldn't be dropped from the geo.simmap to match the data until the LH calculation step (after VCV calculation)
		#o3<-fit_t_generalBETA(drop.tip.simmap(out$geo.simmap,out$geo.simmap$tip.label[which(!out$geo.simmap$tip.label%in%names(M))]), M, fun=new_list_function,diagnostic=T,echo=T)

		DDexpgeo_mvMORPH.lnL<-o3$LogLik
		DDexpgeo_mvMORPH.sig2<-o3$rates[2,1]
		DDexpgeo_mvMORPH.r<-o3$rates[1,1]
		DDexpgeo_mvMORPH.z0<-o3$anc
		DDexpgeo_mvMORPH.conv<-o3$convergence
		DDexpgeo_mvMORPH.time<-proc.time()-pt

		int<-c("fruit",i,names(subdata)[j],length(subtree$tip.label),BM.log_lik,BM.sig2,BM.z0,DDexpgeo.iter,DDexpgeo_RPANDA.lnL,DDexpgeo_RPANDA.sig2,DDexpgeo_RPANDA.r,DDexpgeo_RPANDA.z0,DDexpgeo_RPANDA.conv,DDexpgeo_RPANDA.time[3],DDexpgeo_mvMORPH.lnL,DDexpgeo_mvMORPH.sig2,DDexpgeo_mvMORPH.r,DDexpgeo_mvMORPH.z0,DDexpgeo_mvMORPH.conv,DDexpgeo_mvMORPH.time[3])

		res.mat[counter,]<-int
		print(int)
		write.csv(res.mat,file="tanager_fits_fruit_TESTS.csv")
		counter=counter+1
		}
		print(i)
}
resm<-data.frame(res.mat)
write.csv(resm,file="tanager_fits_fruit_TESTS.csv")
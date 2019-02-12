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

rmats<-resortSMatrix(phylo,r.matrix)
ddm.ob<-createModel_DDexp_multi(phylo,rmats)
params0<-c(0,log(sqrt(var(trait[1,])/max(nodeHeights(phylo)))),-1e-4,-1e-4)	
o1a<-fitTipData(ddm.ob,trait[1,],params0,GLSstyle=TRUE)

#$value
#[1] -42.75546
#
#$inferredParams
#          m0    logsigma0           r1           r2 
#-0.009016524 -1.724855606 -0.046974344 -0.071106515 
#
#$convergence
#[1] 0

ddl.ob<-createModel_DDlin_multi(phylo,rmats)
o1al<-fitTipData(ddl.ob,trait[1,],params0,GLSstyle=TRUE)

#$value
#[1] -42.63477
#
#$inferredParams
#           m0     logsigma0            r1            r2 
#-0.0066771301 -1.8346771805 -0.0006239160 -0.0008672482 
#
#$convergence
#[1] 0


source('~/ddexp/R/fit_t_general_options_old_JC_120219.R', chdir = TRUE)
source('~/ddexp/R/generalized_functions.R')
source('~/Dropbox/Scripts/R scripts/CreateClassObject.R', chdir = TRUE)
source('~/Dropbox/Scripts/R scripts/trimSimmap.R', chdir = TRUE)
source('~/ddexp/R/CreateClassbyClassObject_mvMORPH.R')


phylo=phylo
data=trait[1,]
model="exponential"
regime.map=smap
beta=NULL
sigma=NULL
method=c("Nelder-Mead")
upper=Inf
lower=-Inf ##NOTE this change--this makes all the difference for some reason
control=list(maxit=20000)
diagnostic=TRUE
echo=TRUE
error=NULL

		#first check that regime.map and phylo and data are concordant
		if(!all(as.phylo(phylo)$tip.label == as.phylo(regime.map)$tip.label)) { stop("regime map doesn't match phylogeny")}
		if(length(data) != length(as.phylo(regime.map)$tip.label)) { stop("number of lineages in data and regime map don't match")}
		if(! all (names(data) %in% as.phylo(regime.map)$tip.label)) { stop("names of lineages in data and regime map don't match")}
		if(! all (as.phylo(regime.map)$tip.label %in% names(data)) ) { stop("names of lineages in data and regime map don't match")}
		
		class.object<-try(CreateClassObject(regime.map))
		if(class(class.object)=="try-error"){class.object<-try(CreateClassObject(regime.map,rnd=6))}
		if(class(class.object)=="try-error"){class.object<-CreateClassObject(regime.map,rnd=7)}

		class.df<-return.class.df_subgroup(regime.map,class.object)
		new_list_function<-create.function.list(regime.map,class.object,class.df)
		
		#calculate maxN if DDlin, set to NA if DDexp
		
		maxN<-ifelse(model=="linear",max(class.df[,-1]),NA)
		
		#fit model
		sigma.constraint<-rep(1, dim(regime.map$mapped.edge)[2])
		beta.constraint<-seq(1,by=1,length.out=dim(regime.map$mapped.edge)[2])
		
		o1b<-fit_t_general(tree=regime.map,data=data,fun=new_list_function,error=error, sigma=sigma, beta=beta, model=model,method=method,maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))


#> o1b
#$LogLik
#[1] 42.76091
#
#$AIC
#[1] -77.52181
#
#$AICc
#[1] -76.63292
#
#$rates
#             temp        trop
#beta  -0.07114978 -0.04709504
#sigma  0.03180509  0.03180509
#
#$anc
#[1] -0.009014882
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
#[1] -0.07114978 -0.04709504 -3.44812900
#
#$phyloTrans
#NULL

> -getDataLikelihood(ddm.ob, trait[1,], params=c(o1b$anc,log(sqrt(o1b$rates[2,1])),o1b$rates[1,2],o1b$rates[1,1]))
#[1] 42.75545




		o1bl<-fit_t_general(tree=regime.map,data=data,fun=new_list_function,error=error, sigma=sigma, beta=beta, model="linear",method=method,maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
#$LogLik
#[1] 42.64289
#
#$AIC
#[1] -77.28578
#
#$AICc
#[1] -76.39689
#
#$rates
#               temp          trop
#beta  -0.0008692097 -0.0006261159
#sigma  0.0255324914  0.0255324914
#
#$anc
#[1] -0.00665851
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
#[1] -0.0008692097 -0.0006261159 -3.6678034642
#
#$phyloTrans
#NULL

-getDataLikelihood(ddl.ob, trait[1,], params=c(o1bl$anc,log(sqrt(o1bl$rates[2,1])),o1bl$rates[1,2],o1bl$rates[1,1]))
#[1] 42.63474


###two-slope version is working fine for both cases

######################################
######################################
######################################
######################################

###now back to one slope version

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
source('~/Dropbox/Scripts/R scripts/CreateSMatrix.R')

r.matrix= CreateSMatrix(CreateClassObject(smap),S.cats=c("trop","temp"))
trait<-try(sim_multiDD(smap,pars=matrix(c(0.05,-0.05,-0.1,0),nrow=1),plot=F,S.matrix=r.matrix,rnd=5,model="exponential"))


o1a<-fit_t_comp(phylo,trait[1,],model="DDexp")
#$LH
#[1] 42.41896
#
#$aic
#[1] -78.83791
#
#$aicc
#[1] -78.31617
#
#$free.parameters
#[1] 3
#
#$sig2
#[1] 0.02799225
#
#$r
#[1] -0.02453237
#
#$z0
#[1] -0.006367221
#
#$convergence
#[1] 0

source('~/ddexp/R/fit_t_general_options_old_JC_120219.R', chdir = TRUE)
source('~/ddexp/R/generalized_functions.R')

phylo=phylo
data=trait[1,]
model="exponential"
regime.map=smap
beta=NULL
sigma=NULL
method=c("Nelder-Mead")
upper=Inf
lower=-Inf ##NOTE this change--this makes all the difference for some reason
control=list(maxit=20000)
diagnostic=TRUE
echo=TRUE
error=NULL


		#convert phylo to simmap object
		hold<-rep("A",length(phylo$tip.label))
		hold[1]<-"B"
		names(hold)<-phylo$tip.label
		smap<-make.simmap(phylo,hold,message=F)
		new.maps<-list()
		for(i in 1:length(phylo$edge.length)){
			new.maps[[i]]<-phylo$edge.length[i]
			names(new.maps[[i]])<-"A"
			}
		new.mapped.edge<- as.matrix(rowSums(smap$mapped.edge))
		colnames(new.mapped.edge)<-"A"	
		smap$maps<-new.maps
		smap$mapped.edge<-new.mapped.edge
	
		nodeDist<-vector(mode = "numeric", length = smap$Nnode)
  		totlen<-length(smap$tip.label)
  		root <-totlen  + 1
  		heights<-nodeHeights(smap)
  		for (i in 1:dim(smap$edge)[1]){
    		nodeDist[[smap$edge[i, 1] - totlen]] <- heights[i]
		  }
		nodeDist<-c(nodeDist,max(heights))
		smap$times<-nodeDist
		
		#create function
		class.df<-return.class.df_sympatric(smap)
		new_list_function<-create.function.list_sympatric(smap,class.df)
		
		#calculate maxN if DDlin, set to NA if DDexp
		maxN<-ifelse(model=="linear",max(class.df[,2]),NA)
		
		#fit model
		sigma.constraint<-rep(1, dim(smap$mapped.edge)[2])
		beta.constraint<-rep(1, dim(smap$mapped.edge)[2])

		o1b<-fit_t_general(tree=smap,data=data,fun=new_list_function,error=error, sigma=sigma, beta=beta, model=model,method=method,maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))


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
data<-rnorm(length(BGB.examples$Canidae.phylo$tip.label))
names(data)<-BGB.examples$Canidae.phylo$tip.label

DDexp.geo.fit<-fit_t_comp(BGB.examples$Canidae.phylo, data, model="DDexp", geography.object=Canidae.geography.object)

#$LH
#[1] -58.17102
#
#$aic
#[1] 122.342
#
#$aicc
#[1] 123.142
#
#$free.parameters
#[1] 3
#
#$sig2
#[1] 0.242883
#
#$r
#[1] 0.1209826
#
#$z0
#[1] 0.3174215
#
#$convergence
#[1] 0


DDlin.geo.fit<-fit_t_comp(BGB.examples$Canidae.phylo, data,model="DDlin", geography.object=Canidae.geography.object)

#$LH
#[1] -56.80092
#
#$aic
#[1] 119.6018
#
#$aicc
#[1] 120.4018
#
#$free.parameters
#[1] 3
#
#$sig2
#[1] 2.444631e-07
#
#$b
#[1] 0.08467057
#
#$z0
#[1] 0.3502035
#
#$convergence
#[1] 0


###checking fit_t_DD
source('~/ddexp/R/make.simmap.BGB.R', chdir = TRUE)
source('~/ddexp/R/fit_t_general_options_old_JC_120219.R', chdir = TRUE)
source('~/ddexp/R/generalized_functions.R')


require(geiger)
require(data.table)
require(phytools)
geo.map=make.simmap.BGB(BGB.examples$Canidae.phylo,BGB.examples$Canidae.phylo,ana.events=BGB.examples$Canidae.ana.events,clado.events=BGB.examples$Canidae.clado.events,rnd=6)$geo.simmap

model="exponential"
beta=NULL
sigma=NULL
method=c("Nelder-Mead")
upper=Inf
lower=-Inf ##NOTE this change--this makes all the difference for some reason
control=list(maxit=20000)
diagnostic=TRUE
echo=TRUE
error=NULL


		geo.simmap<-geo.map
		hold<-CreateClassObject(geo.simmap)
		geo.class.df<-return.class.df(geo.simmap,hold)
		new_list_function <- create.function.list(geo.simmap=geo.simmap,geo.class.object=hold, geo.class.df=geo.class.df)

		maxN<-ifelse(model=="linear",max(class.df[,-1]),NA)

		sigma.constraint<-rep(1, dim(geo.map$mapped.edge)[2])
		beta.constraint<-rep(1, dim(geo.map$mapped.edge)[2])

#exponential:
		o1b<-fit_t_general(tree=geo.map,data=data,fun=new_list_function,error=error, sigma=sigma, beta=beta, model=model,method=method,maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))

# successful convergence of the optimizer 
#a reliable solution has been reached 
#
#Summary results for the exponential  model 
#LogLikelihood: 	 -58.17919 
#AIC: 	 122.3584 
#AICc: 	 123.1584 
#3 parameters
#Estimated rates matrix 
#              F        DF         D        BD       BCD      BCDH       BDH     BCDGH      CDGH         C        FI         I
#beta  0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303
#sigma 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937
#            FEI        BI       ABI      ABCI     ABCDI    ABCDFI        CD       ACD      ABCD        CH         B        AB
#beta  0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303
#sigma 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937
#            BFG        BJ        FG         G        BF       BIJ        GH       FGH        FH       CFH      CDFH        AD
#beta  0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303 0.1204303
#sigma 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937 0.2443937
#            CDH       FIJ      AFIJ
#beta  0.1204303 0.1204303 0.1204303
#sigma 0.2443937 0.2443937 0.2443937
#
#Estimated ancestral state 
#0.3172609

#linear:

		o1bl<-fit_t_general(tree=geo.map,data=data,fun=new_list_function,error=error, sigma=sigma, beta=beta, model='linear',method=method,maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))

# successful convergence of the optimizer 
#a reliable solution has been reached 
#
#Summary results for the linear  model 
#LogLikelihood: 	 -56.80715 
#AIC: 	 119.6143 
#AICc: 	 120.4143 
#3 parameters
#Estimated rates matrix 
#                 F           DF            D           BD          BCD         BCDH          BDH        BCDGH         CDGH
#beta  8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02
#sigma 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08
#                 C           FI            I          FEI           BI          ABI         ABCI        ABCDI       ABCDFI
#beta  8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02
#sigma 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08
#                CD          ACD         ABCD           CH            B           AB          BFG           BJ           FG
#beta  8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02
#sigma 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08
#                 G           BF          BIJ           GH          FGH           FH          CFH         CDFH           AD
#beta  8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02 8.471264e-02
#sigma 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08 3.636867e-08
#               CDH          FIJ         AFIJ
#beta  8.471264e-02 8.471264e-02 8.471264e-02
#sigma 3.636867e-08 3.636867e-08 3.636867e-08
#
#Estimated ancestral state 
#0.3510397


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
source('~/Dropbox/Temperate:Tropical/analyses/model_development/DD multi-slope models/DDexpMulti_geo_ADiag.R')
source('~/Dropbox/Temperate:Tropical/analyses/model_development/DD multi-slope models/DDlinMulti_geo_ADiag.R')
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


source('~/ddexp/R/fit_t_general_options_old_JC_120219.R', chdir = TRUE)
source('~/ddexp/R/generalized_functions.R')

phylo=smap
data=M
model="exponential"
subgroup="fruit"
subgroup.map=smap
beta=NULL
sigma=NULL
method=c("Nelder-Mead")
upper=Inf
lower=-Inf ##NOTE this change--this makes all the difference for some reason
control=list(maxit=20000)
diagnostic=TRUE
echo=TRUE
error=NULL

	 	#first check  subgroup map, and phylo and data are concordant
		
		if(!all(phylo$tip.label %in% as.phylo(subgroup.map)$tip.label)) { stop("some lineages in phylogeny don't appear in subgroup map")}
		if( is.null(subgroup) || (!subgroup%in%colnames(subgroup.map$mapped.edge))){ stop("specify a subgroup that appears as a mapped regime in subgroup.map")}

###is this really what we want for subgroup? might inadvertently get rid of trimclass species? DOUBEL CHECK		
		#trim subgroup simmap to species found in focal phylogeny
		subgroup.trimmed<-drop.tip.simmap(subgroup.map,subgroup.map$tip.label[which(!subgroup.map$tip.label%in%phylo$tip.label)])
		
		#trim this diet*region simmap to contain only branches in subgroup (see Drury et al. 2018 PLOS Biol for description of algorithm)
		trimclass.subgroup.trimmed<-trimSimmap(subgroup.trimmed,trim.class=subgroup)

		#create 'class object' for subsequent scripts
		class.object<-try(CreateClassObject(trimclass.subgroup.trimmed))
		
		#add some try catches to prevent crashes due to rounding issues
		if(class(class.object)=="try-error"){class.object<-try(CreateClassObject(trimclass.subgroup.trimmed,rnd=6))}
		if(class(class.object)=="try-error"){class.object<-CreateClassObject(trimclass.subgroup.trimmed,rnd=7)}

		#create class.df that has the diversity of lineages in each simmap class
		subgroup.class.df<-return.class.df_subgroup(trimclass.subgroup.trimmed,class.object)
		
		#set lineages not in subgroup to evolve via BM
		subgroup.class.df[,which(colnames(trimclass.subgroup.trimmed$mapped.edge)!=subgroup)+1]=1 
	
		#trim the subgroup.map to just lineages that are inv at the tips, adjust class.df accordingly
		
		trimclass.subgroup.trimmed.tips<-drop.tip.simmap(trimclass.subgroup.trimmed,trimclass.subgroup.trimmed$tip.label[which(!trimclass.subgroup.trimmed$tip.label%in%names(data))])
		subgroup.class.df.trimmed<-subgroup.class.df[,c(1,match(colnames(trimclass.subgroup.trimmed.tips$mapped.edge),colnames(trimclass.subgroup.trimmed$mapped.edge))+1)]		
	
		subgroup.map.region.root=max(nodeHeights(trimclass.subgroup.trimmed))
		trimclass.subgroup.trimmed.tips.root=max(nodeHeights(trimclass.subgroup.trimmed.tips))

		###	adjust the class.df and class.object if the trimmed tree has a younger root than the ancestral tree
		if(round(subgroup.map.region.root,5)!=round(trimclass.subgroup.trimmed.tips.root,5)){

			trimmed.class.object<-try(CreateClassObject(trimclass.subgroup.trimmed.tips))
			if(class(trimmed.class.object)=="try-error"){trimmed.class.object<-try(CreateClassObject(trimclass.subgroup.trimmed.tips,rnd=6))}
			if(class(trimmed.class.object)=="try-error"){trimmed.class.object<-CreateClassObject(trimclass.subgroup.trimmed.tips,rnd=7)}

			shifted.times<-trimmed.class.object$times+(subgroup.map.region.root-trimclass.subgroup.trimmed.tips.root)
			new.subgroup.class.df.trimmed<-subgroup.class.df.trimmed[c(which(round(out$times,5)==round(min(shifted.times),5)):dim(subgroup.class.df.trimmed)[1]),]
			
			new.subgroup.class.df.trimmed$interval<-c(1:dim(new.subgroup.class.df.trimmed)[1])
			class.object$times<-class.object$times[which(round(class.object$times,5)>=round(min(shifted.times),5))]-round(subgroup.map.region.root-trimclass.subgroup.trimmed.tips.root,5)
			#forces time to start at root of trimmed tree; would be better to pass times directly to new_list_function to avoid overwriting this slot of the out object, which could lead to errors
			subgroup.class.df.trimmed<-new.subgroup.class.df.trimmed
		}
		
		new_list_function<-create.function.list(trimclass.subgroup.trimmed.tips,class.object,subgroup.class.df.trimmed)
		
		#calculate maxN if DDlin, set to NA if DDexp
		
		maxN<-ifelse(model=="linear",max(subgroup.class.df.trimmed[,-1]),NA)
		
		sigma.constraint<-rep(1, dim(trimclass.subgroup.trimmed.tips$mapped.edge)[2])
		beta.constraint<-rep(NA, dim(trimclass.subgroup.trimmed.tips$mapped.edge)[2])
		beta.constraint[which(colnames(trimclass.subgroup.trimmed.tips$mapped.edge)==subgroup)]<-1
		
o3E<-fit_t_general(tree=trimclass.subgroup.trimmed.tips, data=data, fun=new_list_function,error=error, sigma=sigma, beta=beta, model=model,method=method,maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))

# successful convergence of the optimizer 
#a reliable solution has been reached 
#
#Summary results for the exponential  model 
#LogLikelihood: 	 -19.63071 
#AIC: 	 45.26141 
#AICc: 	 45.68999 
#3 parameters
#Estimated rates matrix 
#            fruit        omn        inv
#beta  -0.01412492 0.00000000 0.00000000
#sigma  0.01106355 0.01106355 0.01106355
#
#Estimated ancestral state 
#3.103856

o3L<-fit_t_general(tree=trimclass.subgroup.trimmed.tips, data=data, fun=new_list_function,error=error, sigma=sigma, beta=beta, model="linear",method=method,maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))

# successful convergence of the optimizer 
#a reliable solution has been reached 
#
#Summary results for the linear  model 
#LogLikelihood: 	 -19.77547 
#AIC: 	 45.55094 
#AICc: 	 45.97952 
#3 parameters
#Estimated rates matrix 
#              fruit         omn         inv
#beta  -7.860514e-05 0.000000000 0.000000000
#sigma  9.612849e-03 0.009612849 0.009612849
#
#Estimated ancestral state 
#3.099135


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
source('~/ddexp/R/generalized_functions.R', chdir = TRUE)
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

o1e<-fit_t_comp_subgroup(smap,M,trim.class="fruit",regime.map=smap2,model="DDexp")
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

o1l<-fit_t_comp_subgroup(smap,M,trim.class="fruit",regime.map=smap2,model="DDlin")

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

source('~/ddexp/R/fit_t_general_options_old_JC_120219.R', chdir = TRUE)
source('~/ddexp/R/generalized_functions.R')

phylo=smap
data=M
model="exponential"
subgroup="fruit"
subgroup.map=smap
regime.map=smap2
beta=NULL
sigma=NULL
method=c("Nelder-Mead")
upper=Inf
lower=-Inf ##NOTE this change--this makes all the difference for some reason
control=list(maxit=20000)
diagnostic=TRUE
echo=TRUE
error=NULL


		if(!all(phylo$tip.label %in% as.phylo(subgroup.map)$tip.label)) { stop("some lineages in phylogeny don't appear in subgroup map")}
		if( is.null(subgroup) || (!subgroup%in%colnames(subgroup.map$mapped.edge))){ stop("specify a subgroup that appears as a mapped regime in subgroup.map")}

###This is the version from fits for latitude project	
		#trim subgroup.map and regime.map to species found in focal phylogeny (i.e., a phylogeny trimmed to the taxa of interest)
###is this really what we want for subgroup? might inadvertently get rid of trimclass species? DOUBEL CHECK (I think it's ok for regime map)
		
		subgroup.simmap.trimmed<-drop.tip.simmap(subgroup.map,subgroup.map$tip.label[which(!subgroup.map$tip.label%in%phylo$tip.label)])

		regime.simmap.trimmed<-drop.tip.simmap(regime.map,regime.map$tip.label[which(!regime.map$tip.label%in%phylo$tip.label)])

		##prepare simmap and fit DDexp to subgroup *regime

		#need to build a class.df where diversity represents the number of species in a regime that are ALSO in the subgroup

		class.by.class.object<-try(CreateClassbyClassObject_mvMORPH(map.guild=subgroup.simmap.trimmed,map.regime=regime.map,trim.class=subgroup))
		if(class(class.by.class.object)=="try-error"){class.by.class.object<-try(CreateClassbyClassObject_mvMORPH(map.guild=subgroup.simmap.trimmed,map.regime=regime.map,trim.class=subgroup,rnd=6))}
		if(class(class.by.class.object)=="try-error"){class.by.class.object<-CreateClassbyClassObject_mvMORPH(map.guild=subgroup.simmap.trimmed,map.regime=regime.map,trim.class=subgroup,rnd=7)}
		regime.class.df<-return.class.df_subgroup(class.by.class.object$regime.simmap,class.by.class.object$regime.class.object)

		#set lineages not in "inv" subgroup to evolve via BM
		regime.class.df[,which(colnames(class.by.class.object$regime.simmap$mapped.edge)=='Z')+1]=1

#fix naming scheme to avoid 'region', update annotation accordingly
		#trim the inv region simmap to just lineages that are inv at the tips, adjust class.df accordingly
		regime.simmap.region.trimmed<-drop.tip.simmap(class.by.class.object$regime.simmap,class.by.class.object$regime.simmap$tip.label[which(!class.by.class.object$regime.simmap$tip.label%in%names(data))])
		
		regime.class.df.trimmed<-regime.class.df[,c(1,match(colnames(regime.simmap.region.trimmed$mapped.edge),colnames(class.by.class.object$regime.simmap$mapped.edge))+1)]		

		regime.simmap.region.root=max(nodeHeights(class.by.class.object$regime.simmap))
		regime.simmap.region.trimmed.root=max(nodeHeights(regime.simmap.region.trimmed))

		###	adjust the class.df and class.object if the trimmed tree has a younger root than the ancestral tree
		if(round(regime.simmap.region.root,5)!=round(regime.simmap.region.trimmed.root,5)){
		
			trimmed.class.object<-try(CreateClassObject(regime.simmap.region.trimmed,rnd=5))
			if(class(trimmed.class.object)=="try-error"){trimmed.class.object<-try(CreateClassObject(regime.simmap.region.trimmed,rnd=6))}
			if(class(trimmed.class.object)=="try-error"){trimmed.class.object<-CreateClassObject(regime.simmap.region.trimmed,rnd=7)}

			shifted.times<-trimmed.class.object$times+(regime.simmap.region.root-regime.simmap.region.trimmed.root)
			new.regime.class.df.trimmed<-regime.class.df.trimmed[c(which(round(class.by.class.object$regime.class.object$times,5)==round(min(shifted.times),5)):dim(regime.class.df.trimmed)[1]),]
			
			new.regime.class.df.trimmed$interval<-c(1:dim(new.regime.class.df.trimmed)[1])
			class.by.class.object$regime.class.object$times<-class.by.class.object$regime.class.object$times[which(round(class.by.class.object$regime.class.object$times,5)>=round(min(shifted.times),5))]-round(regime.simmap.region.root-regime.simmap.region.trimmed.root,5)
			#forces time to start at root of trimmed tree; would be better to pass times directly to new_list_function to avoid overwriting this slot of the class.by.class.object object, which could lead to errors
			regime.class.df.trimmed<-new.regime.class.df.trimmed
		}
		
		new_list_function<-create.function.list(regime.simmap.region.trimmed,class.by.class.object$regime.class.object,regime.class.df.trimmed)

		maxN<-ifelse(model=="linear",max(regime.class.df.trimmed[,-1]),NA)
		
		sigma.constraint<-rep(1, dim(regime.simmap.region.trimmed$mapped.edge)[2])
		beta.constraint<-rep(NA, dim(regime.simmap.region.trimmed$mapped.edge)[2])
		beta.constraint[which(colnames(regime.simmap.region.trimmed$mapped.edge)!="Z")]<-1:(dim(regime.simmap.region.trimmed$mapped.edge)[2]-1)
		
	o3e<-fit_t_general(tree=regime.simmap.region.trimmed, data=data, fun=new_list_function,error=error, sigma=sigma, beta=beta, model=model,method=method,maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))

# successful convergence of the optimizer 
#a reliable solution has been reached 
#
#Summary results for the exponential  model 
#LogLikelihood: 	 -19.53093 
#AIC: 	 47.06187 
#AICc: 	 47.78914 
#4 parameters
#Estimated rates matrix 
#        temperate          Z    tropical
#beta  -0.02649129 0.00000000 -0.02673028
#sigma  0.01113273 0.01113273  0.01113273
#
#Estimated ancestral state 
#3.096731

	o3l<-fit_t_general(tree=regime.simmap.region.trimmed, data=data, fun=new_list_function,error=error, sigma=sigma, beta=beta, model="linear",method=method,maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))


# successful convergence of the optimizer 
#a reliable solution has been reached 
#
#Summary results for the linear  model 
#LogLikelihood: 	 -20.46343 
#AIC: 	 48.92686 
#AICc: 	 49.65413 
#4 parameters
#Estimated rates matrix 
#         temperate           Z     tropical
#beta  2.834454e-05 0.000000000 2.297838e-07
#sigma 5.891906e-03 0.005891906 5.891906e-03
#
#Estimated ancestral state 
#3.089635


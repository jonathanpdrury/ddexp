##testing subgroup pruning algorithm

load('~/ddexp/data/tanager.data.Rd')

require(deSolve)
require(geiger)
require(phytools)
library(corpcor)
require(methods)
require(tools)
require(RPANDA)
require(data.table)

source('~/ddexp/R/fit_t_general_subgroup.R')
source('~/ddexp/R/generalized_functions.R')
source('~/ddexp/R/CreateGeobyClassObject_mvMORPH.R')
source('~/ddexp/R/make.simmap.BGB.R')
source('~/PANDA/R/fit_t_comp_subgroup.R')
source('~/Dropbox/Scripts/R scripts/trimSimmap.R', chdir = TRUE)

tan.tree<-tanager.data[[1]]
tandata<-read.csv("~/Dropbox/Manuscripts/PLOS Biology Tanager Project/data/MASTER_tandata.csv")
 
guild.map<-tanager.data[[4]]


fruit<-subset(tandata, diet=="fruit")#62 spp, #59 in tree
mass<-fruit[,5]
names(mass)<-fruit[,3]
nc<-name.check(tan.tree,mass)
fruit.tree<-drop.tip(tan.tree,nc$tree_not_data)
subdata<-fruit[which(fruit$treename%in%fruit.tree$tip.label),]


j=6
i=1

#first trim data and subtree to frugivores, then fit BM model

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

#then fit subgroup DDexp model in RPANDA
		
iter=i
guild.mapi<-drop.tip.simmap(guild.map[[i]],guild.map[[i]]$tip.label[which(!guild.map[[i]]$tip.label%in%tan.tree$tip.label)])
pt<-proc.time()
o5<-try(fit_t_comp_subgroup(full.phylo=tan.tree,ana.events=tanager.data[[2]][[1]][[i]],clado.events=tanager.data[[2]][[2]][[i]],map=guild.mapi,data=M,trim.class="fruit",model="DDexp",par=NULL,method="Nelder-Mead",bounds=NULL))
while(class(o5)=="try-error"){
	iter=iter+1
	guild.mapi<-drop.tip.simmap(guild.map[[iter]],guild.map[[iter]]$tip.label[which(!guild.map[[iter]]$tip.label%in%tan.tree$tip.label)])
	o5<-try(fit_t_comp_subgroup(full.phylo=tan.tree,ana.events=tanager.data[[2]][[1]][[i]],clado.events=tanager.data[[2]][[2]][[i]],map=guild.mapi,data=M,trim.class="fruit",model="DDexp",par=NULL,method="Nelder-Mead",bounds=NULL))
	}
DDexpgeo.iter<-iter	
DDexpgeo_RPANDA.lnL<-o5$LH
DDexpgeo_RPANDA.sig2<-o5$sig2
DDexpgeo_RPANDA.r<-o5$r
DDexpgeo_RPANDA.z0<-o5$z0
DDexpgeo_RPANDA.conv<-o5$convergence

DDexpgeo_RPANDA.time<-proc.time()-pt

###now fit subgroup DDexp model in mvMORPH

out<-CreateGeobyClassObject_mvMORPH(phylo=tan.tree,map=guild.mapi,ana.events=tanager.data[[2]][[1]][[i]],clado.events=tanager.data[[2]][[2]][[i]],trim.class="fruit")
geo.class.df<-return.class.df(out$geo.simmap,out$geo.class.object)
geo.class.df[,which(colnames(out$geo.simmap$mapped.edge)=='Z')+1]=0
##to check against BM, can use this approach: 		geo.class.df[,2:dim(geo.class.df)[2]][!is.na(geo.class.df[,2:dim(geo.class.df)[2]])]=0

new_list_function<-create.function.list(out$geo.simmap,out$geo.class.object,geo.class.df)

full.simmap.tree=out$geo.simmap
data=M
fun=new_list_function
error=NULL 
beta=NULL
sigma=NULL
model=c("exponential")
method=c("L-BFGS-B")
upper=Inf
lower=-20
control=list(maxit=20000)
diagnostic=TRUE
echo=TRUE
constraint=TRUE

#o3<-try(fit_t_general_subgroup(full.simmap.tree=out$geo.simmap,data= M, fun=new_list_function,diagnostic=T,echo=T))
##NOTE: this is where the problem was, tips shouldn't be dropped from the geo.simmap to match the data until the LH calculation step (after VCV calculation)
#o3<-fit_t_generalBETA(drop.tip.simmap(out$geo.simmap,out$geo.simmap$tip.label[which(!out$geo.simmap$tip.label%in%names(M))]), M, fun=new_list_function,diagnostic=T,echo=T)

#tree=drop.tip.simmap(out$geo.simmap,out$geo.simmap$tip.label[which(!out$geo.simmap$tip.label%in%names(M))])
#data=M
#fun=new_list_function
#error=NULL 
#beta=NULL
#sigma=NULL
#model=c("exponential")
#method=c("L-BFGS-B")
#upper=Inf
#lower=-20
#control=list(maxit=20000)
#diagnostic=TRUE
#echo=TRUE
#constraint=TRUE
#
#DDexpgeo_mvMORPH.lnL<-o3$LogLik
#DDexpgeo_mvMORPH.sig2<-o3$rates[2,1]
#DDexpgeo_mvMORPH.r<-o3$rates[1,1]
#DDexpgeo_mvMORPH.z0<-o3$anc
#DDexpgeo_mvMORPH.conv<-o3$convergence
#DDexpgeo_mvMORPH.time<-proc.time()-pt
#
#int<-c("fruit",i,names(subdata)[j],length(subtree$tip.label),BM.log_lik,BM.sig2,BM.z0,DDexpgeo.iter,DDexpgeo_RPANDA.lnL,DDexpgeo_RPANDA.sig2,DDexpgeo_RPANDA.r,DDexpgeo_RPANDA.z0,DDexpgeo_RPANDA.conv,DDexpgeo_RPANDA.time[3],DDexpgeo_mvMORPH.lnL,DDexpgeo_mvMORPH.sig2,DDexpgeo_mvMORPH.r,DDexpgeo_mvMORPH.z0,DDexpgeo_mvMORPH.conv,DDexpgeo_mvMORPH.time[3])
#
#print(int)

#I now think that the problem with this wasn't dropping the tips, per se, but dropping the tips without adjusting the geo.class.df to be ordered correctly. Let's see if this fixes things:

geo.simmap.trimmed<-drop.tip.simmap(out$geo.simmap,out$geo.simmap$tip.label[which(!out$geo.simmap$tip.label%in%names(M))])
geo.class.df.trimmed<-geo.class.df[,c(1,match(colnames(geo.simmap.trimmed$mapped.edge),colnames(out$geo.simmap$mapped.edge))+1)]
#class.object.trimmed<-CreateClassObject(geo.simmap.trimmed)


##ISSUE with this approach (passing modified class.df to function), I think, is that the time values for looking up states with geo.class.df no longer works if the trimmed phylogeny root is not also the root of the full phylogeny


geo.simmap.root=max(nodeHeights(out$geo.simmap))
geo.simmap.trimmed.root=max(nodeHeights(geo.simmap.trimmed))

if(round(geo.simmap.root,5)!=round(geo.simmap.trimmed.root,5)){stop("root of trimmed tree is younger than root of trimclass tree")}


new_list_function<-create.function.list(geo.simmap.trimmed,out$geo.class.object,geo.class.df.trimmed)
source('~/ddexp/R/fit_t_general.R')

o3<-fit_t_general(geo.simmap.trimmed, M, fun=new_list_function,diagnostic=T,echo=T)



###running tanager tests on my computer:

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
		o5<-try(fit_t_comp_subgroup(full.phylo=tan.tree,ana.events=tanager.data[[2]][[1]][[i]],clado.events=tanager.data[[2]][[2]][[i]],map=guild.mapi,data=M,trim.class="fruit",model="DDexp",par=NULL,method="Nelder-Mead",bounds=NULL))
		while(class(o5)=="try-error"){
			iter=iter+1
			guild.mapi<-drop.tip.simmap(guild.map[[iter]],guild.map[[iter]]$tip.label[which(!guild.map[[iter]]$tip.label%in%tan.tree$tip.label)])
			o5<-try(fit_t_comp_subgroup(full.phylo=tan.tree,ana.events=tanager.data[[2]][[1]][[i]],clado.events=tanager.data[[2]][[2]][[i]],map=guild.mapi,data=M,trim.class="fruit",model="DDexp",par=NULL,method="Nelder-Mead",bounds=NULL))
			}
		DDexpgeo.iter<-iter	
		DDexpgeo_RPANDA.lnL<-o5$LH
		DDexpgeo_RPANDA.sig2<-o5$sig2
		DDexpgeo_RPANDA.r<-o5$r
		DDexpgeo_RPANDA.z0<-o5$z0
		DDexpgeo_RPANDA.conv<-o5$convergence
		
		DDexpgeo_RPANDA.time<-proc.time()-pt
		
		###now add in fit_t_* tools and record output as well as time
		out<-CreateGeobyClassObject_mvMORPH(phylo=tan.tree,map=guild.mapi,ana.events=tanager.data[[2]][[1]][[i]],clado.events=tanager.data[[2]][[2]][[i]],trim.class="fruit")
		geo.class.df<-return.class.df(out$geo.simmap,out$geo.class.object)
		#geo.class.df[,which(colnames(out$geo.simmap$mapped.edge)=='Z')+1]=0
		geo.class.df[,which(colnames(out$geo.simmap$mapped.edge)=='Z')+1]=1 #this is the change from the 0802 tests
		##to check against BM, can use this approach: 		geo.class.df[,2:dim(geo.class.df)[2]][!is.na(geo.class.df[,2:dim(geo.class.df)[2]])]=0

		geo.simmap.trimmed<-drop.tip.simmap(out$geo.simmap,out$geo.simmap$tip.label[which(!out$geo.simmap$tip.label%in%names(M))])
		geo.class.df.trimmed<-geo.class.df[,c(1,match(colnames(geo.simmap.trimmed$mapped.edge),colnames(out$geo.simmap$mapped.edge))+1)]
		#class.object.trimmed<-CreateClassObject(geo.simmap.trimmed)
		
		
		
		geo.simmap.root=max(nodeHeights(out$geo.simmap))
		geo.simmap.trimmed.root=max(nodeHeights(geo.simmap.trimmed))
		
		if(round(geo.simmap.root,5)!=round(geo.simmap.trimmed.root,5)){stop("root of trimmed tree is younger than root of trimclass tree")}
		
		new_list_function<-create.function.list(geo.simmap.trimmed,out$geo.class.object,geo.class.df.trimmed)
		

		pt<-proc.time()
		#o3<-fit_t_general_subgroup(out$geo.simmap, M, fun=new_list_function,diagnostic=T,echo=T)
		#o3<-fit_t_generalBETA(drop.tip.simmap(out$geo.simmap,out$geo.simmap$tip.label[which(!out$geo.simmap$tip.label%in%names(M))]), M, fun=new_list_function,diagnostic=T,echo=T)
		o3<-fit_t_general(geo.simmap.trimmed, M, fun=new_list_function,diagnostic=T,echo=T)
		
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
#write.csv(resm,file="tanager_fits_fruit_TESTS_20180802.csv")
write.csv(resm,file="tanager_fits_fruit_TESTS_20180804.csv")


##that looks a *lot* better; question is still whether this is a robust solution; that is, whether root of trimmed.simmap is likely to be older than the root of the clade--maybe ignore for now
##biggest difference is with i=1, j=2 (0.060898463 = delta lnL). Could this be because of additional trimming owing to missing data?

##TURNS OUT this is because I was turning allopatric lineages into n=0 rather than n=1; now that I've fixed this, the delta lnLs between the two approaches max out at 0.009605526

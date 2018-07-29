
require(RPANDA)
require(phytools)
load('anole.data.Rd')

##first create bank of trees
simmap.list<-list()

state.list<-list()
state.list[[1]]<-c(rep("A",50),rep("B",50))
state.list[[2]]<-c(rep("A",25),rep("B",25),rep("C",25),rep("D",25))
state.list[[3]]<-c(rep("A",10),rep("B",10),rep("C",10),rep("D",10),rep("E",10),rep("F",10),rep("G",10),rep("H",10),rep("I",10),rep("J",10))

counter=1
for(i in 1:length(state.list)){
	#create simmap with ordered, mutually exclusive ranges
	states<-state.list[[i]]
	names(states)<-anole.data[[1]]$tip.label
	simmap1<-make.simmap(anole.data[[1]],states,model="ARD")
	simmap.list[[counter]]<-simmap1
	counter=counter+1
	
	#create simmap with randomized, mutually exclusive ranges
	states<-sample(state.list[[i]])
	names(states)<-anole.data[[1]]$tip.label
	simmap2<-make.simmap(anole.data[[1]],states,model="ARD")
	simmap.list[[counter]]<-simmap2
	counter=counter+1

	if(i >1){
	#create simmap with ordered, nonmutually exclusive ranges
	
	states<-state.list[[i]]
	states[which(states=="B")]<-"AC"
	names(states)<-anole.data[[1]]$tip.label
	simmap3<-make.simmap(anole.data[[1]],states,model="ARD")
	simmap.list[[counter]]<-simmap3
	counter=counter+1

	#create simmap with randomized, nonmutually exclusive ranges
	states<-sample(state.list[[i]])
	states[which(states=="B")]<-"AC"
	names(states)<-anole.data[[1]]$tip.label
	simmap4<-make.simmap(anole.data[[1]],states,model="ARD")
	simmap.list[[counter]]<-simmap4
	counter=counter+1

	if(i==3){
	states<-state.list[[i]]
	states[which(states=="B")]<-"AC"
	states[which(states=="E")]<-"DF"
	states[which(states=="H")]<-"GI"

	names(states)<-anole.data[[1]]$tip.label
	simmap5<-make.simmap(anole.data[[1]],states,model="ARD")
	simmap.list[[counter]]<-simmap5
	counter=counter+1
	
	states<-sample(state.list[[i]])
	states[which(states=="B")]<-"AC"
	states[which(states=="E")]<-"DF"
	states[which(states=="H")]<-"GI"

	names(states)<-anole.data[[1]]$tip.label
	simmap6<-make.simmap(anole.data[[1]],states,model="ARD")
	simmap.list[[counter]]<-simmap6
	counter=counter+1
	}
	}
}

save(simmap.list,file="simmap.list.RData")

res.mat<-matrix(nrow=length(simmap.list),ncol=16)
colnames(res.mat)<-c("i","no.states","mutually.exclusive","dim.class.df","lnL.RPANDA","sig2.RPANDA","slope.RPANDA","z0.RPANDA","convergence.RPANDA","time.RPANDA","lnL.mvMORPH","sig2.mvMORPH","slope.mvMORPH","z0.mvMORPH","convergence.mvMORPH","time.mvMORPH")
mexclusive<-c(1,1,1,1,0,0,1,1,0,0,0,0)
source('~/Dropbox/Scripts/R scripts/CreateGeoObject_overlap_v2.R', chdir = TRUE)

for(i in 11:length(simmap.list)){
	if(i %in% c(1,3,5,7,9,11:12)){
	anole.geo.object<-CreateGeoObject_overlap(anole.data[[1]],simmap.list[[i]])
	pc<-proc.time()
	o1<-fit_t_comp(anole.data[[1]],anole.data[[3]],model="DDexp",geography.object=anole.geo.object)
	time.RPANDA<-proc.time()-pc
	lnL.RPANDA<-o1$LH
	sig2.RPANDA<-o1$sig2
	slope.RPANDA<-o1$r
	z0.RPANDA<-o1$z0
	convergence.RPANDA<-o1$convergence
	
	
	class.object<-CreateClassObject(simmap.list[[i]],rnd=7)
	geo.class.df<-return.class.df(simmap.list[[i]],class.object)
	new_list_function <- create.function.list(geo.simmap=simmap.list[[i]],geo.class.object=class.object, class.df=geo.class.df)
	
	pc<-proc.time()
	o2<-fit_t_general(simmap.list[[i]],anole.data[[3]],fun=new_list_function)
	time.mvMORPH<-proc.time()-pc
	lnL.mvMORPH<-o2$LogLik
	sig2.mvMORPH<-o2$rates[2,1]
	slope.mvMORPH<-o2$rates[1,1]
	z0.mvMORPH<-o2$anc
	convergence.mvMORPH<-o2$convergence
	int<-c(i,dim(geo.class.df)[2],mexclusive[i],dim(geo.class.df)[1],lnL.RPANDA,sig2.RPANDA,slope.RPANDA,z0.RPANDA,convergence.RPANDA,time.RPANDA[3],lnL.mvMORPH,sig2.mvMORPH,slope.mvMORPH,z0.mvMORPH,convergence.mvMORPH,time.mvMORPH[3])
	print(int)
	res.mat[i,]<-int
	}	
}

##all fits that I was able to run had almost identical likelihoods, including non-mutually exclusive and rapidly evolving states; 
##so discrepancy between approaches for tanagers doesn't appear to be an issue with the way event timings happen or implementation of shared range states
##other possibilities include:
#something about subgroup trimming? the simmap production might be off?
#something about issues arising when there are 56 class combinations? 
#could try to implement ranges as mutually exclusive to see if that helps
#timing is way faster in this approach too, yet it is slower in tanager stuff



##TANAGERS testing with complex biogeography, 100 tips culled from tree
load('tanager.data.Rd')


tanager.tree.trimmed<-drop.tip(tanager.data[[1]],tanager.data[[1]]$tip.label[101:323])
tan.geo.object<-CreateBioGeoB_Object_subclade(anc.phylo=tanager.data[[1]],subclade.phylo= tanager.tree.trimmed,ana.events=tanager.data[[2]][[1]][[1]],clado.events=tanager.data[[2]][[2]][[1]])
tanmass.trimmed<-tanager.data[[3]][which(names(tanager.data[[3]])%in%tanager.tree.trimmed$tip.label)]

out<-make.simmap.BGB(anc.phylo=tanager.data[[1]],subclade.phylo=tanager.tree.trimmed,ana.events=tanager.data[[2]][[1]][[1]],clado.events=tanager.data[[2]][[2]][[1]])
geo.simmap<-out$geo.simmap
class.object<-out$class.object

geo.class.df<-return.class.df(geo.simmap,class.object)
new_list_function <- create.function.list(geo.simmap=out$geo.simmap,geo.class.object=class.object, class.df=geo.class.df)

pc<-proc.time()
fit_t_general(out$geo.simmap,tanmass.trimmed,fun=new_list_function)
proc.time()-pc


##try using geo.object tools for fit_t_comp
tan.geo.object<-CreateGeoObject_overlap(tanager.tree.trimmed,geo.simmap)


pc<-proc.time()
o1<-fit_t_comp(tanager.tree,tanmass.trimmed,model="DDexp",geography.object=tan.geo.object)
proc.time()-pc

##I was hoping to build a geo.object using the same simmap passed to fit_t_general--if it returned the same likelihood, this would suggest an issue with make.simmap.BGB; but
##it looks like the geo.object isn't working, so I need to revisit CreateGeoObject_overlap_v2.R;
##maybe add a pass through the reduction step? but that doesn't explain so many repeats of time intervals



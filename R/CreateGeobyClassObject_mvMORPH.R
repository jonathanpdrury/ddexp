#This is a wrapper script that merges a class object and a geographic stochastic map from BioGeoBEARS 
##OUTPUTS geo.simmap, geo.class.object, and subgroup.class.object
##in geo.class.object, range is replaced with 'allo' if it is not in subgroup identified by trim class at that time step

#The inputs are:
#phylo (phylogeny used to build ancestral range stochastic map) 
#ana.events (anagenetic events tables from ancestral range stochastic map); if stratified BGB output, must be output from stratified_BGB_to_tables.R
#clado.events (cladogenetic events tables from ancestral range stochastic map); if stratified BGB output, must be output from stratified_BGB_to_tables.R
#map (simmap object that is a stochastic map of discrete traits used to identify guild of competitors) 
#trim.class (discrete category identifying species in 'data' (e.g., herbivores))

CreateGeobyClassObject_mvMORPH<-function(phylo,map,trim.class,ana.events,clado.events,rnd=5){

trc=trim.class

##trim tree

new.map<-trimSimmap(map,trc)

##create class.object
class.object<-CreateClassObject(new.map)

##update geo.object
#geo.object<-CreateBioGeoB_Object_subclade(phylo,new.map,ana.events,clado.events)
##need this to return geo.simmap, too

out<-make.simmap.BGB(phylo,new.map,ana.events,clado.events)
geo.simmap<-out$geo.simmap
geo.class.object<-out$class.object

##first concatenate geo.object and class.object timings 
	geot<-round(geo.class.object$times,rnd)
	clat<-round(class.object$times,rnd)
	nodeDist<-sort(unique(c(geot,clat)))
	nodeDiff<-diff(nodeDist)
	if(any(nodeDiff<= (2*(10^-rnd)))){stop("potential rounding error, two time bins very similar, try changing rnd digits")}
	

#initialize counter for geo.class.object and class.object
u<-0 #class object
y<-0 #geo object



###NEED TO UPDATE THIS so that only class objects are reconciled, not multiplied against one another
hold.class<-list()
hold.geo<-list()

for(i in 1:length(nodeDiff)){

	if((nodeDist[i]%in%geot) && (nodeDist[i]%in%clat)){ #if timing is the same for both
		u = u+1
		y = y+1
		class.int<-class.object$class.object[[u]]
		geo.int<-geo.class.object$class.object[[y]]
		geo.int[which(class.int[,2]!=trc),2]<-'allo'
		hold.class[[i]]<-class.int
		hold.geo[[i]]<-geo.int
	}
	if((nodeDist[i]%in%geot) && (!nodeDist[i]%in%clat)){ #this means that geo.object changes but class object doesn't
		y = y+1
		class.int<-class.object$class.object[[u]]
		geo.int<-geo.class.object$class.object[[y]]
		geo.int[which(class.int[,2]!=trc),2]<-'allo'
		hold.class[[i]]<-class.int
		hold.geo[[i]]<-geo.int
	}
	if((!nodeDist[i]%in%geot) && (nodeDist[i]%in%clat)){ #this means that class.object changes but geo object doesn't
		u = u+1
		class.int<-class.object$class.object[[u]]
		geo.int<-geo.class.object$class.object[[y]]
		geo.int[which(class.int[,2]!=trc),2]<-'allo'
		hold.class[[i]]<-class.int
		hold.geo[[i]]<-geo.int
	}
	}

	
return(list(subgroup.simmap=new.map,subgroup.class.object=list(class.object=hold.class,times=nodeDist,spans=nodeDiff),geo.simmap=geo.simmap,geo.class.object=list(class.object=hold.geo,times=nodeDist,spans=nodeDiff)))#new phylo object, #new times, #new spans, #new geo object
}
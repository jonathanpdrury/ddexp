#This is a wrapper script that merges a class object and a geographic stochastic map from BioGeoBEARS 
##OUTPUTS geo.simmap, geo.class.object, and subgroup.class.object
##in geo.class.object, range is replaced with 'Z' if it is not in subgroup identified by trim class at that time step

#The inputs are:
#phylo (phylogeny used to build ancestral range stochastic map) 
#ana.events (anagenetic events tables from ancestral range stochastic map); if stratified BGB output, must be output from stratified_BGB_to_tables.R
#clado.events (cladogenetic events tables from ancestral range stochastic map); if stratified BGB output, must be output from stratified_BGB_to_tables.R
#map (simmap object that is a stochastic map of discrete traits used to identify guild of competitors) 
#trim.class (discrete category identifying species in 'data' (e.g., herbivores))

CreateGeobyClassObject_mvMORPH<-function(phylo,map,trim.class,ana.events,clado.events,trimmed.class.object=NULL,trimmed.geo.simmap.output=NULL,rnd=5){

trc=trim.class

##trim tree

new.map<-trimSimmap(map,trc)

##create class.object if not provided

if(is.null(trimmed.class.object)){
class.object<-CreateClassObject(new.map,rnd=rnd)
}else{
class.object<-trimmed.class.object
}

if(is.null(trimmed.geo.simmap.output)){
out<-make.simmap.BGB(phylo,new.map,ana.events,clado.events,return.mat=TRUE)
} else{
out<-trimmed.geo.simmap.output
if(is.null(out$mat)){stop("'mat' slot required from geo.simmap")}
if(is.null(out$geo.simmap)){stop("'geo.simmap' slot required from geo.simmap")}
if(is.null(out$class.object)){stop("'class.object' slot required from geo.simmap")}
}

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
		geo.int[which(class.int[,2]!=trc),2]<-'Z'
		hold.class[[i]]<-class.int
		hold.geo[[i]]<-geo.int
	}
	if((nodeDist[i]%in%geot) && (!nodeDist[i]%in%clat)){ #this means that geo.object changes but class object doesn't
		y = y+1
		class.int<-class.object$class.object[[u]]
		geo.int<-geo.class.object$class.object[[y]]
		geo.int[which(class.int[,2]!=trc),2]<-'Z'
		hold.class[[i]]<-class.int
		hold.geo[[i]]<-geo.int
	}
	if((!nodeDist[i]%in%geot) && (nodeDist[i]%in%clat)){ #this means that class.object changes but geo object doesn't
		u = u+1
		class.int<-class.object$class.object[[u]]
		geo.int<-geo.class.object$class.object[[y]]
		geo.int[which(class.int[,2]!=trc),2]<-'Z'
		hold.class[[i]]<-class.int
		hold.geo[[i]]<-geo.int
	}
	}


##update geo.simmap to add 'Z' to allopatric times


	phylo<-geo.simmap
	mat<-out$mat
	maps.list=list()
	
	for(k in 1:length(phylo$edge.length)){
		
		#identify branch from edge matrix
		#lookup the name of this branch in the 'mat' matrix compiled above
		#lookup which nat elements have the name of this branch
		#write a vector of the nodeDiff values named with the ranges for each of these elements

		
		lf<-phylo$edge[k,1]
		ri<-phylo$edge[k,2]
		br<-mat[which(mat[,1]==lf & mat[,3]==ri),2]
		natis<-which(sapply(hold.geo,function(x)br%in%x[,1]))
		out.vec<-nodeDiff[natis]
		name.vec<-vector()
		for(n in 1:length(out.vec)){
			name.vec<-c(name.vec,hold.geo[[natis[n]]][which(hold.geo[[natis[n]]][,1]==br),2])	
		}	
		names(out.vec)<-name.vec
		
		#sum adjacent elements with the same name
		
		out.vec.simple<-vector()
		counter=1
		for(i in 1: length(out.vec)){
			if(i == 1 || i == (counter+1)){
				while((length(out.vec)>counter) && (names(out.vec[i])==names(out.vec[counter+1]))){
					counter=counter+1
					}
				hold<-sum(out.vec[i:counter])
				names(hold)<-names(out.vec[i])
				out.vec.simple<-c(out.vec.simple,hold)	
				}
		}
		
		
		maps.list[[k]]<-out.vec.simple
	
	}

	mapped.edge<-matrix(nrow=dim(phylo$edge)[1],ncol=length(unique(names(unlist(maps.list)))))
	colnames(mapped.edge)<-unique(names(unlist(maps.list)))
	
	for(k in 1:dim(phylo$edge)[1]){
		for(j in 1:dim(mapped.edge)[2]){
			hold<-which(names(maps.list[[k]])==colnames(mapped.edge)[j])
			mapped.edge[k,j]<-ifelse(length(hold)==0,0,sum(maps.list[[k]][hold]))
		}
	}
	
	outsmap<-list(edge=phylo$edge,edge.length=phylo$edge.length,tip.label=phylo$tip.label,Nnode=phylo$Nnode,maps=maps.list,mapped.edge=mapped.edge,Q="NA",logL="NA")
	class(outsmap)<-c("phylo","simmap")

	
return(list(subgroup.simmap=new.map,subgroup.class.object=list(class.object=hold.class,times=nodeDist,spans=nodeDiff),geo.simmap=outsmap,geo.class.object=list(class.object=hold.geo,times=nodeDist,spans=nodeDiff)))#new phylo object, #new times, #new spans, #new geo object
}
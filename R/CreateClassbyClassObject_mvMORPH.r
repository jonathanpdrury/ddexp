#This is a wrapper script that merges two class objects
##OUTPUTS geo.simmap, geo.class.object, and subgroup.class.object
##in geo.class.object, range is replaced with 'Z' if it is not in subgroup identified by trim class at that time step

#The inputs are:
#map.guild (simmap object that is a stochastic map of discrete traits used to identify guild of competitors) (e.g., diet category)
#map.regime (simmap object that is a stochastic map of discrete traits used to identify the regime of interest) (e.g., tropical vs. temperate)
#trim.class (discrete category identifying species in 'data' (e.g., herbivores))

CreateClassbyClassObject_mvMORPH<-function(map.guild,map.regime,trim.class,rnd=5){

trc=trim.class

##trim tree

new.map<-trimSimmap(map.guild,trc)

##create class.object if not provided

class.object.guild<-CreateClassObject(new.map,rnd=rnd)

#trim regime to match guild

map.regime.trimmed<-drop.tip.simmap(map.regime,map.regime$tip.label[which(!map.regime$tip.label%in%new.map$tip.label)])

class.object.regime<-CreateClassObject(map.regime.trimmed,rnd=rnd,return.mat=TRUE)

##first concatenate geo.object and class.object timings 
	regt<-round(class.object.regime$times,rnd)
	guit<-round(class.object.guild$times,rnd)
	nodeDist<-sort(unique(c(regt,guit)))
	nodeDiff<-diff(nodeDist)
	if(any(nodeDiff<= (2*(10^-rnd)))){stop("potential rounding error, two time bins very similar, try changing rnd digits")}
	

#initialize counter for geo.class.object and class.object
u<-0 #guild class object
y<-0 #regime class object


###NEED TO UPDATE THIS so that only class objects are reconciled, not multiplied against one another
hold.guild<-list()
hold.regime<-list()

for(i in 1:length(nodeDiff)){

	if((nodeDist[i]%in%regt) && (nodeDist[i]%in%guit)){ #if timing is the same for both
		u = u+1
		y = y+1
		gui.int<-class.object.guild$class.object[[u]]
		reg.int<-class.object.regime$class.object[[y]]
		reg.int[which(gui.int[,2]!=trc),2]<-'Z'
		hold.guild[[i]]<-gui.int
		hold.regime[[i]]<-reg.int
	}
	if((nodeDist[i]%in%regt) && (!nodeDist[i]%in%guit)){ #this means that geo.object changes but class object doesn't
		y = y+1
		gui.int<-class.object.guild$class.object[[u]]
		reg.int<-class.object.regime$class.object[[y]]
		reg.int[which(gui.int[,2]!=trc),2]<-'Z'
		hold.guild[[i]]<-gui.int
		hold.regime[[i]]<-reg.int
	}
	if((!nodeDist[i]%in%regt) && (nodeDist[i]%in%guit)){ #this means that class.object changes but geo object doesn't
		u = u+1
		gui.int<-class.object.guild$class.object[[u]]
		reg.int<-class.object.regime$class.object[[y]]
		reg.int[which(gui.int[,2]!=trc),2]<-'Z'
		hold.guild[[i]]<-gui.int
		hold.regime[[i]]<-reg.int
	}
	}


	phylo<-class.object.regime$phylo
	mat<-class.object.regime$mat
	maps.list=list()
	
	for(k in 1:length(phylo$edge.length)){
		
		#identify branch from edge matrix
		#lookup the name of this branch in the 'mat' matrix compiled above
		#lookup which nat elements have the name of this branch
		#write a vector of the nodeDiff values named with the ranges for each of these elements

		
		lf<-phylo$edge[k,1]
		ri<-phylo$edge[k,2]
		br<-mat[which(mat[,1]==lf & mat[,3]==ri),2]
		natis<-which(sapply(hold.regime,function(x)br%in%x[,1]))
		out.vec<-nodeDiff[natis]
		name.vec<-vector()
		for(n in 1:length(out.vec)){
			name.vec<-c(name.vec,hold.regime[[natis[n]]][which(hold.regime[[natis[n]]][,1]==br),2])	
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

	
return(list(subgroup.simmap=new.map,subgroup.class.object=list(class.object=hold.guild,times=nodeDist,spans=nodeDiff),regime.simmap=outsmap,regime.class.object=list(class.object=hold.regime,times=nodeDist,spans=nodeDiff)))#new phylo object, #new times, #new spans, #new geo object

}
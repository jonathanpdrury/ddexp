#Version September 06,2016

#' Simulates data under a model of matching competition with two competitive regimes, specified in S.matrix
#'
#' @param phylo A phylogenetic tree
#' @param pars A matrix, columns containing values for sig2, S1 (first row of S.matrix), S2 (second row of S.matrix), S3 (third row of S.matrix), root.value
#' @param OU OU or BM model
#' @param min.Nsegments Min number of segments for the simulation. Default is 2500
#' @param plot should the last run be plotted?
#' @param S.matrix 
#' @param allo.index; parameter number of the S category that indicates 'allopatric' for subgroup trimming (sets diversity to 1)

#' @return A matrix with the
#' @details sensitivity to interspecific similarity is defined by alpha parameter maximum amount of repulsion is defined by the max parameter
#' @export
#' @examples #need to update examples with biogeo.object

sim_multiDD_subgroup<-function(phylo,pars, verbose=TRUE, min.Nsegments=2500, plot=FALSE,S.matrix,allo.index=4,rnd=6){
  paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
  paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS
	#something to check that it isn't sorted
  nodeDist<-vector(mode = "numeric", length = phylo$Nnode)
  root <- length(phylo$tip.label) + 1
  heights<-phytools::nodeHeights(phylo)
  totlen<-max(heights)
  len<-root-1
  nodeDist<-c(as.numeric(sort(max(ape::branching.times(phylo))-ape::branching.times(phylo))),totlen)
  nodeDiff<-diff(nodeDist)
  if(!is.null(phylo$node.label)){phylo$node.label<-NULL}
  old.labels<-as.numeric(names(sort(ape::branching.times(phylo),decreasing=TRUE)))
  old.edge<-phylo$edge
  if(any(diff(old.labels)!=1)){ #if nodes are not in sequential order, this renames them so that they are
  	checkmat<-cbind(old.labels,seq(root,len+phylo$Nnode))
  	old.edge<-phylo$edge
  	for(j in 1:phylo$Nnode){phylo$edge[which(old.edge==checkmat[j,1])]<-checkmat[j,2]}
  	}
  mat<-matrix(nrow=0, ncol=3)
  counter_three_letters <- 0
  for(i in 1:phylo$Nnode){
    other<-phylo$edge[phylo$edge[,1]==i+len, 2]
    for(b in other){
      int<-matrix(ncol=3)
      int[1]<-i+len
      if(b>len){
        counter_three_letters <- counter_three_letters + 1
        int[2]<-paste(".",THREELETTERS[counter_three_letters],sep="")
        int[3]<-b
      } else {
        int[2]<-phylo$tip.label[[b]]
        int[3]<-0 ##NOTE :I am considering tips to be "0" here and use this below
      }
      mat<-rbind(mat,int)
    }
  }
  nat<-list()
  branchesPresent <- rep(NA, length(nodeDiff))
  for(i in 1:length(nodeDiff)){
    if(i==1){
        nat[[i]]<-mat[mat[,1]==(len+i),2]
      } else {
        IN<-vector()
        P<-mat[as.numeric(mat[,1])<=(len+i),c(2,3)]
        IN<-c(IN, P[P[,2]=="0",1],P[as.numeric(P[,2])>(len+i),1])
        nat[[i]]<-IN
      }
    branchesPresent[i] = length(nat[[i]])
  }
	
  masterbranch<-list()
  segmentSize <- rep(NA, phylo$Nnode)
  mappings <- list()
  segsize = sum(nodeDiff)/min.Nsegments
  if(segsize>min(S.matrix$spans)){print(paste0('segsize=',segsize,'; consider increasing N segments'))}

  for(i in 1:phylo$Nnode){ ##for each node interval

    nati<-nat[[i]]
    mappings[[i]] = rep(NA, branchesPresent[i])

    if(i==1){
      segmentSize[i] = ceiling(nodeDiff[i] / segsize)
    }else{
      segmentSize[i] = ceiling(nodeDiff[i] / segsize)
      for (j in 1:branchesPresent[i]){
        if(length(which(nati[j] == pastNati)) > 0){
          mappings[[i]][j] = which(nati[j] == pastNati)
        }else{
          mappings[[i]][j] = which(mat[mat[,3]==mat[mat[,2]==nati[j],1],2] == pastNati)
        }
      }
    }
    pastNati <- nati
  }

  ## Looping over the parameters

  out <- matrix( nrow = nrow(pars), ncol = len)

  for (p in 1:nrow(pars)){

    sig2 = pars[p,1]
    S1 = pars[p,2]
    S2 = pars[p,3]
    
    if(allo.index!=4){stop("re-create Smatrix, setting the third row to allopatric species for subgroup trimming")}
    S3 = pars[p,4]
    
    root.value = pars[p,5]

	smat = S.matrix$S.matrix
	if(!all(colnames(smat[[length(smat)]])==nat[[length(nat)]])){stop("order of S.matrix species is incorrect; perhaps this is a sorted matrix?")}	
    newDist<-S.matrix$times
    newDiff<-S.matrix$spans

	if((sig2/(sig2*exp(S1*length(phylo$tip.label))))>=1e3){warning("r1 parameter leads to sig2 at tips that is 1000x smaller than root value, consider changing")}
	if((sig2/(sig2*exp(S2*length(phylo$tip.label))))>=1e3){warning("r2 parameter leads to sig2 at tips that is 1000x smaller than root value, consider changing")}


    timecount=1
    for(i in 1:phylo$Nnode){
        traitMat <- matrix(nrow = branchesPresent[i], ncol = segmentSize[i]+1)

        if (i == 1){
          traitMat[,1] = root.value
        }else{
          traitMat[,1] = masterbranch[[i-1]][mappings[[i]],(segmentSize[i-1]+1)] #added +1 here since segmentSize[i-1] is penultimate column, not the last one
        }

        tempInd<- 1:branchesPresent[i] # hack to have fast selection of not k, seemed to be faster than a call to which()

        for(k in 1:segmentSize[i]){

	          	for(j in 1:branchesPresent[i]){
	            
	            if(dim(smat[[timecount]])[2]!=branchesPresent[i]){stop(paste0("error on i =",i," k=",k," j=",j))}
       		 	tempS1 =  ifelse(smat[[timecount]][1,j]>0,exp(S1*sum(smat[[timecount]][1,])),0)
       		 	tempS2 =  ifelse(smat[[timecount]][2,j]>0,exp(S2*sum(smat[[timecount]][2,])),0)
       		 	tempS3 =  ifelse(smat[[timecount]][3,j]>0,exp(S3*1),0)
       		 	
            	traitMat[j,k+1]<-traitMat[j,k] +rnorm(1,0,sqrt(sig2*((tempS1+tempS2+tempS3)/sum(smat[[timecount]][,j]))*segsize))

				if((k!=segmentSize[i])&&(j==branchesPresent[i])){timecount= ifelse(round(((k*segsize)+nodeDist[i]),rnd)>=round(newDist[timecount+1],rnd),timecount+1,timecount)}
				###loop for last segment size (to preserve exact branch lengths)
				if(k==segmentSize[i] && nodeDiff[i]%%segsize!=0){
					segsizeT= nodeDiff[i]%%segsize
	            	traitMat[j,k+1]<-traitMat[j,k] +rnorm(1,0,sqrt(sig2*((tempS1+tempS2+tempS3)/sum(smat[[timecount]][,j]))*segsizeT))
				
				 	if(j==branchesPresent[i]){timecount= ifelse(round((((k-1)*segsize)+nodeDist[i]+segsizeT),rnd)>=round(newDist[timecount+1],rnd),timecount+1,timecount)}
					}
	        	}	        	
			}
        masterbranch[[i]] = traitMat
      }
    out[p,] = as.vector(masterbranch[[phylo$Nnode]][, tail(segmentSize,n=1)+1])
    if(p==1){colnames(out)<-unlist(nat[[i]])}
  }

  if(plot==TRUE){
    print("plotting last simulated dataset")
    M=seq(0,sum(nodeDiff),length=sum(segmentSize))
    O=list()
    for(i in 1:length(nodeDiff)){
    O[[i]]<-data.frame(seq(nodeDist[[i]],nodeDist[[i+1]],length=(segmentSize[[i]]+1)),as.data.frame(t(masterbranch[[i]])))
    }
    t.plot<-plot(M,1:length(M),col="white", ylim=c(range(sapply(masterbranch,range))), xlab="time", ylab="Value")
    for(i in 1:length(nodeDiff)){
    	for(j in 1:(ncol(O[[i]])-1)){
    		lines(O[[i]][,1],O[[i]][,j+1])
    	}
    }
  }

  return(out)
}



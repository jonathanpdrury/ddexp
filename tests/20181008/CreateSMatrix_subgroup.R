CreateSMatrix_subgroup<-function(class.obj,S.cats=list()){
	if(length(S.cats)!=3){stop("S.cats must contain a list of 3 vectors specifying which states correspond to which slope regimes in the DDM model")}
	S1=S.cats[[1]]
	S2=S.cats[[2]]
	S3=S.cats[[3]]
	out.list<-list()
	for(i in 1:length(class.obj$class.object)){
		imat<-class.obj$class.object[[i]]
		out.mat<-matrix(nrow=3,ncol=dim(imat)[1])
		colnames(out.mat)<-imat[,1]
		out.mat[1,]<-imat[,2]%in%S1
		out.mat[2,]<-imat[,2]%in%S2
		out.mat[3,]<-imat[,2]%in%S3
		out.list[[i]]<-out.mat*1
		}
	return(list(S.matrix=out.list,times=class.obj$times,spans=class.obj$spans))
	}	
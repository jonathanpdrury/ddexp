CreateSMatrix<-function(class.obj,S.cats=list()){
	if(length(S.cats)!=2){stop("S.cats must contain a list of two vectors specifying which states correspond to which S values in the multiS matching competition model")}
	S1=S.cats[[1]]
	S2=S.cats[[2]]
	out.list<-list()
	for(i in 1:length(class.obj$class.object)){
		imat<-class.obj$class.object[[i]]
		out.mat<-matrix(nrow=2,ncol=dim(imat)[1])
		colnames(out.mat)<-imat[,1]
		out.mat[1,]<-imat[,2]%in%S1
		out.mat[2,]<-imat[,2]%in%S2
		out.list[[i]]<-out.mat*1
		}
	return(list(S.matrix=out.list,times=class.obj$times,spans=class.obj$spans))
	}	
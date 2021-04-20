return_target<-function(tree,target,start.slice,slice.interval){
	slice=start.slice
	hold<-treeSlice(tree,slice=slice)
	while(any(sapply(hold,function(hold)(length(hold$tip.label)))>target)){
		slice=slice+slice.interval
		hold<-treeSlice(tree,slice=slice)
		}
print(slice)	
return(hold)
}	
	
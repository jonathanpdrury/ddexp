##tests on 20170731; updated make.simmap.BGB to have simplified maps (adjacent elements with same name are summed)

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

#> fit_t_general(out$geo.simmap,tanmass.trimmed,fun=new_list_function)
#
# successful convergence of the optimizer 
#a reliable solution has been reached 
#
#Summary results for the exponential  model 
#LogLikelihood: 	 -69.56709 
#AIC: 	 145.1342 
#AICc: 	 145.3842 
#3 parameters
#Estimated rates matrix 
#             [,1]
#beta  -0.09160398
#sigma  0.19332660
#
#Estimated ancestral state 
#2.81236
#> proc.time()-pc
#   user  system elapsed 
#  6.947   0.095   7.117 


##try using geo.object tools for fit_t_comp
tan.geo.object<-CreateGeoObject_overlap(tanager.tree.trimmed,geo.simmap)


pc<-proc.time()
o1<-fit_t_comp(tanager.tree.trimmed,tanmass.trimmed,model="DDexp",geography.object=tan.geo.object)
proc.time()-pc

#   user  system elapsed 
# 58.986   0.994  60.219 
#> 
#> 
#> o1
#$LH
#[1] -71.50201
#
#$aic
#[1] 149.004
#
#$aicc
#[1] 149.254
#
#$free.parameters
#[1] 3
#
#$sig2
#[1] 0.2308094
#
#$r
#[1] -0.08876915
#
#$z0
#[1] 2.814415
#
#$convergence
#[1] 0

##Using this approach of building simmap and then geo.object from simmap, I get the exactly same likelihood
##as I get when I build the geo.object directly from BGB objects; shows that make.simmap.BGB and CreateGeoObject_overlap
##are NOT the source of the discrepancy between the two approaches

#which leaves the following possibilities:
#(1) an error with scripts generating class.df or functions to pass to fit_t_general()
##a general error might be caught by treating ranges as mutually exclusive and comparing likelihoods
#if there is still a discrepancy, then error is with one of these functions
#if not, then error could be in how class.df is defined for non mutually exclusive taxa

#(2) an error with something in fit_t_general related to non-mutual exclusivity (but what could this be??)


##testing possibility 1 

return.class.df_mutually.exclusive<-function(simmap,class.object){
	states<-colnames(simmap$mapped.edge)
	for(i in 1:length(states)){
		eval(parse(text=paste('d',i,'<-sapply(class.object$class.object,function(x)sum(x[,2]==states[',i,']))',sep="")))
	}
	eval(parse(text=paste('return(data.frame(interval=1:length(d1),',paste('d',1:length(states),sep="",collapse=','),'))',sep="")))
}

tanager.tree.trimmed<-drop.tip(tanager.data[[1]],tanager.data[[1]]$tip.label[101:323])
tan.geo.object<-CreateBioGeoB_Object_subclade(anc.phylo=tanager.data[[1]],subclade.phylo= tanager.tree.trimmed,ana.events=tanager.data[[2]][[1]][[1]],clado.events=tanager.data[[2]][[2]][[1]])
tanmass.trimmed<-tanager.data[[3]][which(names(tanager.data[[3]])%in%tanager.tree.trimmed$tip.label)]

out<-make.simmap.BGB(anc.phylo=tanager.data[[1]],subclade.phylo=tanager.tree.trimmed,ana.events=tanager.data[[2]][[1]][[1]],clado.events=tanager.data[[2]][[2]][[1]])
geo.simmap<-out$geo.simmap
class.object<-out$class.object

geo.class.df<-return.class.df_mutually.exclusive(geo.simmap,class.object)
new_list_function <- create.function.list(geo.simmap=out$geo.simmap,geo.class.object=class.object, class.df=geo.class.df)

pc<-proc.time()
fit_t_general(out$geo.simmap,tanmass.trimmed,fun=new_list_function)
proc.time()-pc
#
# successful convergence of the optimizer 
#a reliable solution has been reached 
#
#Summary results for the exponential  model 
#LogLikelihood: 	 -54.96317 
#AIC: 	 115.9263 
#AICc: 	 116.1763 
#3 parameters
#Estimated rates matrix 
#            [,1]
#beta  -0.3779437
#sigma  0.0657500
#
#Estimated ancestral state 
#2.719351
#> proc.time()-pc
#   user  system elapsed 
#  2.830   0.063   2.892 

tan.geo.object<-RPANDA::CreateGeoObject(tanager.tree.trimmed,geo.simmap)
pc<-proc.time()
fit_t_comp(tanager.tree.trimmed,tanmass.trimmed,model="DDexp",geography.object=tan.geo.object)
proc.time()-pc

#> fit_t_comp(tanager.tree.trimmed,tanmass.trimmed,model="DDexp",geography.object=tan.geo.object)
#$LH
#[1] -57.67175
#
#$aic
#[1] 121.3435
#
#$aicc
#[1] 121.5935
#
#$free.parameters
#[1] 3
#
#$sig2
#[1] 0.0707564
#
#$r
#[1] -0.3102167
#
#$z0
#[1] 2.725467
#
#$convergence
#[1] 0
#
#> proc.time()-pc
#   user  system elapsed 
# 44.615   1.863  46.446 


##Looks like the discrepancy hasn't got anything to do with mutual exclusivity per se, but rather accounting for states
##Could be because:
##(1)ordering of states gets messed up somewhere and so script calls the wrong column

##note that parameter values are very similar between approaches; maybe increasing number of states impacts likelihood calculation somehow?

##different length(class.object$times) and tan.geo.object$times, why?
##looks like CreateClassObject(geo.simmap) returns something different from what comes out of make.simmap.BGB(); suggests that there is an error in make.simmap.BGB when trimming by subclade happens that affects class.object, rather than geo.simmap

out<-make.simmap.BGB(anc.phylo=tanager.data[[1]],subclade.phylo=tanager.tree.trimmed,ana.events=tanager.data[[2]][[1]][[1]],clado.events=tanager.data[[2]][[2]][[1]])
geo.simmap<-out$geo.simmap
class.object<-out$class.object
hold<-CreateClassObject(geo.simmap)
geo.class.df<-return.class.df_mutually.exclusive(geo.simmap,hold)
new_list_function <- create.function.list(geo.simmap=out$geo.simmap,geo.class.object=hold, class.df=geo.class.df)

pc<-proc.time()
fit_t_general(out$geo.simmap,tanmass.trimmed,fun=new_list_function)
proc.time()-pc


# successful convergence of the optimizer 
#a reliable solution has been reached 
#
#Summary results for the exponential  model 
#LogLikelihood: 	 -57.68923 
#AIC: 	 121.3785 
#AICc: 	 121.6285 
#3 parameters
#Estimated rates matrix 
#             [,1]
#beta  -0.30975043
#sigma  0.07070193
#
#Estimated ancestral state 
#2.725821
#> proc.time()-pc
#   user  system elapsed 
#  3.967   0.066   4.027 


##well, that's more like it. trying now with non-mutually exclusive stuff:

geo.class.df<-return.class.df(geo.simmap,hold)
new_list_function <- create.function.list(geo.simmap=out$geo.simmap,geo.class.object=hold, class.df=geo.class.df)

pc<-proc.time()
fit_t_general(out$geo.simmap,tanmass.trimmed,fun=new_list_function)
proc.time()-pc

#> fit_t_general(out$geo.simmap,tanmass.trimmed,fun=new_list_function)
#
# successful convergence of the optimizer 
#a reliable solution has been reached 
#
#Summary results for the exponential  model 
#LogLikelihood: 	 -71.5022 
#AIC: 	 149.0044 
#AICc: 	 149.2544 
#3 parameters
#Estimated rates matrix 
#             [,1]
#beta  -0.08870122
#sigma  0.23048672
#
#Estimated ancestral state 
#2.814249
#> proc.time()-pc
#   user  system elapsed 
#  8.300   0.130   8.304 
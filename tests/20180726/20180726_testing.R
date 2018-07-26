##Double checking that general functions for creating diversity functions scripts work

require(RPANDA)
load('anole.data.Rd')

anole.geo.object<-CreateGeoObject(anole.data[[1]],anole.data[[2]][[1]])

class.object<-CreateClassObject(anole.data[[2]][[1]],rnd=5)
geo.class.df<-return.class.df(anole.data[[2]][[1]],class.object)
new_list_function <- create.function.list(geo.simmap=anole.data[[2]][[1]],geo.class.object=class.object, class.df=geo.class.df)

class(anole.data[[2]][[1]])<-c("phylo",'simmap')
fit_t_general(anole.data[[2]][[1]],anole.data[[3]],fun=new_list_function)

# successful convergence of the optimizer 
#a reliable solution has been reached 
#
#Summary results for the exponential  model 
#LogLikelihood: 	 -8.066213 
#AIC: 	 22.13243 
#AICc: 	 22.38243 
#3 parameters
#Estimated rates matrix 
#              [,1]
#beta  -0.042092236
#sigma  0.008689931
#
#Estimated ancestral state 
#0.001343236


fit_t_comp(anole.data[[1]],anole.data[[3]],model="DDexp",geography.object=anole.geo.object)

#$LH
#[1] -8.073313
#
#$aic
#[1] 22.14663
#
#$aicc
#[1] 22.39663
#
#$free.parameters
#[1] 3
#
#$sig2
#[1] 0.008694668
#
#$r
#[1] -0.04209809
#
#$z0
#[1] 0.0008694877
#
#$convergence
#[1] 0
#
#attr(,"class")
#[1] "fit_t.comp"


###MOTMOT (non-mutually exclusive biogeography) testing

smap<-stratified_BGB_to_tables(motmot.data[[1]],motmot.data[[2]],1)
geo.object<-CreateBioGeoB_Object_subclade(anc.phylo=motmot.data[[1]],subclade.phylo=motmot.data[[1]],ana.events=smap$ana.int,clado.events=smap$clado.int)

out<-make.simmap.BGB(anc.phylo=motmot.data[[1]],subclade.phylo=motmot.data[[1]],ana.events=smap$ana.int,clado.events=smap$clado.int)
geo.simmap<-out$geo.simmap
class.object<-out$class.object

geo.class.df<-return.class.df(geo.simmap,class.object)
new_list_function <- create.function.list(geo.simmap=out$geo.simmap,geo.class.object=class.object, class.df=geo.class.df)

fit_t_general(out$geo.simmap,motmot.data[[3]],fun=new_list_function)

# successful convergence of the optimizer 
#a reliable solution has been reached 
#
#Summary results for the exponential  model 
#LogLikelihood: 	 -6.76888 
#AIC: 	 19.53776 
#AICc: 	 23.53776 
#3 parameters
#Estimated rates matrix 
#             [,1]
#beta  0.098886415
#sigma 0.005700228
#
#Estimated ancestral state 
#4.315413


#> fit_t_comp(motmot.data[[1]],motmot.data[[3]],model="DDexp",geography.object=geo.object)
#$LH
#[1] -6.768879
#
#$aic
#[1] 19.53776
#
#$aicc
#[1] 23.53776
#
#$free.parameters
#[1] 3
#
#$sig2
#[1] 0.005703389
#
#$r
#[1] 0.09879885
#
#$z0
#[1] 4.315408
#
#$convergence
#[1] 0



##testing with tanagers, 100 tips culled from tree


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

# convergence of the optimizer has not been reached, try simpler model 
#a reliable solution has been reached 
#
#Summary results for the exponential  model 
#LogLikelihood: 	 -69.54682 
#AIC: 	 145.0936 
#AICc: 	 145.3436 
#3 parameters
#Estimated rates matrix 
#             [,1]
#beta  -0.09174758
#sigma  0.19397255
#
#Estimated ancestral state 
#2.812528


#> proc.time()-pc
#   user  system elapsed 
#379.976   3.938 388.902 


pc<-proc.time()
fit_t_comp(tanager.tree.trimmed, tanmass.trimmed,model="DDexp",geography.object=tan.geo.object)
proc.time()-pc

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

#> proc.time()-pc
#   user  system elapsed 
# 63.852   0.633  65.558 



###fit_t_general approach takes much, much longer than fit_t_comp approach in this case, and gives the wrong likelihood for that parameter set

likelihood_t_DD_geog(tanager.tree.trimmed,tanmass.trimmed,par=c(log(0.2308094),-0.08876915),tan.geo.object,model="DDexp")
#[1] 71.50201
#^ correct lnL is returned for fit_t_comp

likelihood_t_DD_geog(tanager.tree.trimmed,tanmass.trimmed,par=c(log(0.19397255),-0.09174758),tan.geo.object,model="DDexp")
#[1] 73.21124
#^ incorrect lnL is returned for fit_t_general


##trying script profiling

profFile<-"prof.out"
Rprof(profFile, line.profiling=TRUE,memory.profiling=TRUE)
fit_t_general(out$geo.simmap,tanmass.trimmed,fun=new_list_function)

Rprof(NULL)

profFileout<-summaryRprof(profFile, lines="show", memory="both")
save(profFileout,file="profFileout.RData")


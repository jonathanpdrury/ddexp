#This script has a reproducible example of the optimisation issues with fit_t_general.R


require(phytools)
require(RPANDA)

set.seed(1234)

fiftytipyuletrees<-pbtree(n=50,nsim=1)
regions<-c("trop","trop","trop","temp","temp")
phylo<-fiftytipyuletrees
samp.reg<-sample(regions)
trait<-c(rep(samp.reg[1],10),rep(samp.reg[2],10),rep(samp.reg[3],10),rep(samp.reg[4],10),rep(samp.reg[5],10))
names(trait)<-phylo$tip.label
smap<-make.simmap(phylo,trait,model="ARD")

source('sim_DDmulti.R')
source('DDexpMulti_nogeo_ADiag.R')
source('DDlinMulti_nogeo_ADiag.R')
source('resortSMatrix.R')
source('CreateSMatrix.R')
source('fit_t_general_options_old.R')
source('generalized_functions.R')

r.matrix= CreateSMatrix(CreateClassObject(smap),S.cats=c("trop","temp"))
trait<-try(sim_multiDD(smap,pars=matrix(c(0.05,-0.05,-0.1,0),nrow=1),plot=F,S.matrix=r.matrix,rnd=5,model="exponential"))


###now fitting model with RPANDA tools:

rmats<-resortSMatrix(phylo,r.matrix)
ddl.ob<-createModel_DDlin_multi(phylo,rmats)
params0<-c(0,log(sqrt(var(trait[1,])/max(nodeHeights(phylo)))),-1e-4,-1e-4)	
o1al<-fitTipData(ddl.ob,trait[1,],params0,GLSstyle=TRUE)
o1al


##now fitting model with fit_t_general tools

phylo=phylo
data=trait[1,]
model="linear"
regime.map=smap
beta=NULL
sigma=NULL
method=c("L-BFGS-B")
upper=Inf
lower=-20
control=list(maxit=20000)
diagnostic=TRUE
echo=TRUE
error=NULL

		#first check that regime.map and phylo and data are concordant
		if(!all(as.phylo(phylo)$tip.label == as.phylo(regime.map)$tip.label)) { stop("regime map doesn't match phylogeny")}
		if(length(data) != length(as.phylo(regime.map)$tip.label)) { stop("number of lineages in data and regime map don't match")}
		if(! all (names(data) %in% as.phylo(regime.map)$tip.label)) { stop("names of lineages in data and regime map don't match")}
		if(! all (as.phylo(regime.map)$tip.label %in% names(data)) ) { stop("names of lineages in data and regime map don't match")}
		
		class.object<-try(CreateClassObject(regime.map))
		if(class(class.object)=="try-error"){class.object<-try(CreateClassObject(regime.map,rnd=6))}
		if(class(class.object)=="try-error"){class.object<-CreateClassObject(regime.map,rnd=7)}

		class.df<-return.class.df_subgroup(regime.map,class.object)
		new_list_function<-create.function.list(regime.map,class.object,class.df)
		
		#calculate maxN if DDlin, set to NA if DDexp
		
		maxN<-ifelse(model=="linear",max(class.df[,-1]),NA)
		
		#fit model
		sigma.constraint<-rep(1, dim(regime.map$mapped.edge)[2])
		beta.constraint<-seq(1,by=1,length.out=dim(regime.map$mapped.edge)[2])
		
		out<-fit_t_general(tree=regime.map,data=data,fun=new_list_function,error=error, sigma=sigma, beta=beta, model=model,method=method,maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		print(out)
		
## JC TEST WITH RPANDA FUNCTION
		
		-getDataLikelihood(ddl.ob, trait[1,], params=c(out$anc,log(sqrt(out$rates[2,1])),out$rates[1,2],out$rates[1,1]))
		# ok here it return the same log-likelihood value as fit_t_general. expected since the beta are 0 and the approximations come from here
		
##note that model doesn't converge (convergence message == 52) using  L-BFGS-B or Nelder-Mead

method="BB"

		out<-fit_t_general(tree=regime.map,data=data,fun=new_list_function,error=error, sigma=sigma, beta=beta, model=model,method=method,maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
		print(out)
		
##now model converges fine

##however, the likelihood of this fit is 39.33939, which is close to, but far enough from 39.99859 (from the RPANDA approach above) to be worrying

##looking at the likelihood of the RPANDA version using the ML parameters from the fit_t_general approach:

-getDataLikelihood(ddl.ob, trait[1,], params=c(out$anc,log(sqrt(out$rates[2,1])),out$rates[1,2],out$rates[1,1]))

##returns 39.3478, which is very near the ML estimate from fit_t_general, which leads me to conclude that while the two scripts produce similar likelihoods, the optimiser isn't hitting the ML value in the fit_t_general script

		
		
#### TESTS & MODIFS JC ------------------------------------------------------------------------------


# Fit with BFGS and Nelder-Mead by setting the lower bound to -Inf and removing the part on l. 210-214
# test=sigma+(beta*maxN)
# if(any(test<=0)){
#	LL<-list()
#	LL$logl<-Inf
#	}else{
#
# I also added on line 13 (just in case):
# tree <- reorderSimmap(tree, order="postorder")

# Nelder-Mead
outT1<-fit_t_general(tree=regime.map, data=data,fun=new_list_function, error=error, 
                   sigma=NULL, beta=NULL, model=model, method="Nelder-Mead", maxN=maxN, upper=Inf, 
                   lower=-Inf, control=control, diagnostic=diagnostic, echo=echo, constraint=list(sigma=sigma.constraint, beta=beta.constraint))
print(outT1) # log lik is [1] 43.07889

# BB (longer but provide similar result)
outT2<-fit_t_general(tree=regime.map, data=data,fun=new_list_function, error=error, 
                     sigma=NULL, beta=NULL, model=model, method="BB", maxN=maxN, upper=Inf, 
                     lower=-20, control=control, diagnostic=diagnostic, echo=echo, constraint=list(sigma=sigma.constraint, beta=beta.constraint))
print(outT2)

# L-BFGS-B
outT3<-fit_t_general(tree=regime.map, data=data,fun=new_list_function, error=error, 
                     sigma=NULL, beta=NULL, model=model, method="L-BFGS-B", maxN=maxN, upper=Inf, 
                     lower=-Inf, control=control, diagnostic=diagnostic, echo=echo, constraint=list(sigma=sigma.constraint, beta=beta.constraint))
print(outT3)

# try different starting values?
outT3<-fit_t_general(tree=regime.map, data=data,fun=new_list_function, error=error, 
                     sigma=exp(outT1$param[3]), beta=c(outT1$param[1:2]), model=model, method="L-BFGS-B", maxN=maxN, upper=Inf, 
                     lower=-Inf, control=control, diagnostic=diagnostic, echo=echo, constraint=list(sigma=sigma.constraint, beta=beta.constraint))
print(outT3)


sigma <- outT1$param[3]
beta <- outT1$param[1:2]
outfix<-fit_t_general(tree=regime.map, data=data, fun=new_list_function, error=error, sigma=sigma, beta=beta, model=model,
                    method="fixed", maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
print(outfix)

#### TESTS & MODIFS JC -- fixed parameters instead of estimation to see if I can recover the log-likelihood

sigma <- exp(c(o1al$inferredParams[2]))^2
beta <- o1al$inferredParams[4:3] # REVERT THE PARAM, IT's NOT THE SAME ORDER
out1<-fit_t_general(tree=regime.map, data=data, fun=new_list_function, error=error, sigma=sigma, beta=beta, model=model,
                    method="fixed", maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
print(out1)
# ok very similar results here for the LL	


sigma <- exp(outT1$param[3])
beta <- outT1$param[1:2]
out2<-fit_t_general(tree=regime.map, data=data, fun=new_list_function, error=error, sigma=sigma, beta=beta, model=model,
                    method="fixed", maxN=maxN, upper=upper, lower=lower, control=control,diagnostic=diagnostic, echo=echo,constraint=list(sigma=sigma.constraint, beta=beta.constraint))
print(out2) # just to check


## COMPARE WITH EXP model
outT1e<-fit_t_general(tree=regime.map, data=data,fun=new_list_function, error=error, 
                     sigma=NULL, beta=NULL, model="exponential", method="Nelder-Mead", maxN=maxN, upper=Inf, 
                     lower=-Inf, control=control, diagnostic=diagnostic, echo=echo, constraint=list(sigma=sigma.constraint, beta=beta.constraint))
print(outT1e) # log lik is [1] 44.097
# not that far from the linear model...


# COMPARE to PANDA?
-getDataLikelihood(ddl.ob, trait[1,], params=c(outT1$anc,log(sqrt(exp(outT1$param[3]))),outT1$param[2],outT1$param[1]))
#??? don't know why...
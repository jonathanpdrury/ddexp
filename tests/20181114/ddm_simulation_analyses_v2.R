require(data.table)
require(ggplot2)

simdat<-read.csv("~/ddexp/tests/20181114/ddm_simulation_study_20181114.csv")

par(mfrow=c(2,1))
plot(simdat$sim.sig2_1,simdat$est.sig2_1)
abline(a=0,b=1,col="red")
plot(simdat$sim.sig2_2,simdat$est.sig2_2)
abline(a=0,b=1,col="red")
plot(simdat$sim.sig2_3,simdat$est.sig2_3)
abline(a=0,b=1,col="red")

#need to double check if sig2 is correct across scripts/estimates (i.e., not a log/exp or sqrt transformation?)
##need to check that regime_3 and Z are correct (i.e., that non-group membership was coded right)
##double check simulator
##could parameter values have been set too low?

boxplot(est.r1~factor(sim.r1),data=simdat)
abline(h=c(-0.1,-0.05,-0.025,0,0.05),lty=2,col="red")
boxplot(est.r2~factor(sim.r2),data=simdat)
abline(h=c(-0.1,-0.05,-0.025,0,0.05),lty=2,col="red")
boxplot(est.r1-est.r2~factor(sim.r1-sim.r2),data=simdat,xlab="slope_1 - slope_2 (simulated)", ylab="slope_1 - slope_2 (estimated)")
abline(h=c(-0.1,-0.05,-0.025,0,0.025,0.05,0.1),lty=2,col="red")

simdat<-subset(simdat,remove==0)
dt<-data.table(simdat)
dt[,list(median.ddm.wi=median(DDM.wi),max.ddm.wi=max(DDM.wi)),by="sim.r1,sim.r2"]
models.col=c("#ffb3b3","#cc0000","#00b38f","#009933","#FFD732","#E18700")
names(models.col)<-c("BM1.wi","BMM.wi","OU1.wi","OUM.wi","DD1.wi","DDM.wi")

##come up with pipeline for analyzing model selection
m.dat.1<-melt(simdat,id.vars=c("X","sim.r1","sim.r2"),measure.vars=c("BM1.wi","BMM.wi","OU1.wi","OUM.wi","DD1.wi","DDM.wi"))
ggplot(m.dat.1,aes(x=interaction(sim.r1,sim.r2),y=value,fill=variable))+geom_bar(position="fill",stat="identity")+scale_fill_manual(values=models.col)+theme_bw()


##these results look much, much better, so it appears this model fitting framework is better/more reliable than the unconstrained DDM that estimates a sig2 for each regime, including the catch-all "out of regime" state Z

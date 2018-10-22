require(data.table)
require(ggplot2)

simdat<-read.csv("~/ddexp/tests/20181008/ddm_simulation_study_201810.csv")

plot(simdat$sim.sig2_1,simdat$est.sig2_1)
plot(simdat$sim.sig2_2,simdat$est.sig2_2)
plot(simdat$sim.sig2_3,simdat$est.sig2_3)

#need to double check if sig2 is correct across scripts/estimates (i.e., not a log/exp or sqrt transformation?)
##need to check that regime_3 and Z are correct (i.e., that non-group membership was coded right)
##double check simulator
##could parameter values have been set too low?

boxplot(est.r1~factor(sim.r1),data=simdat)
abline(h=c(-0.1,-0.05,-0.025,0,0.05),lty=2,col="red")
boxplot(est.r2~factor(sim.r2),data=simdat)
abline(h=c(-0.1,-0.05,-0.025,0,0.05),lty=2,col="red")
hist(simdat$est.r3)

#perhaps this is a case for forcing to 0?

dt<-data.table(simdat)
dt[,list(median.ddm.wi=median(DDM.wi),max.ddm.wi=max(DDM.wi)),by="sim.r1,sim.r2"]

##come up with pipeline for analyzing model selection
m.dat.1<-melt(simdat,id.vars=c("X","sim.r1","sim.r2"),measure.vars=c("BM1.wi","BMM.wi","OU1.wi","OUM.wi","DD1.wi","DDM.wi"))
ggplot(m.dat.1,aes(x=interaction(sim.r1,sim.r2),y=value,fill=variable))+geom_bar(position="fill",stat="identity")+theme_bw()
##check which regime has more taxa, could asymmetry in number shed any light on these patterns?

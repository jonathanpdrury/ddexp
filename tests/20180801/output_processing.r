tf<-read.csv("~/ddexp/tests/20180801/tanager_fits_fruit_TESTS.csv")

par(mfrow=c(2,2))
plot(tf$DDexpgeo_RPANDA.lnL,tf$DDexpgeo_mvMORPH.lnL)
abline(a=0,b=1,col="red")
plot(tf$DDexpgeo_RPANDA.r,tf$DDexpgeo_mvMORPH.r)
abline(a=0,b=1,col="red")
plot(sqrt(tf$DDexpgeo_RPANDA.sig2),sqrt(tf$DDexpgeo_mvMORPH.sig2))
abline(a=0,b=1,col="red")
plot(tf$DDexpgeo_RPANDA.z0,tf$DDexpgeo_mvMORPH.z0)
abline(a=0,b=1,col="red")

tf<-read.csv("~/ddexp/tests/20180802/tanager_fits_fruit_TESTS_20180802.csv")

par(mfrow=c(2,2))
plot(tf$DDexpgeo_RPANDA.lnL,tf$DDexpgeo_mvMORPH.lnL)
abline(a=0,b=1,col="red")
plot(tf$DDexpgeo_RPANDA.r,tf$DDexpgeo_mvMORPH.r)
abline(a=0,b=1,col="red")
plot(sqrt(tf$DDexpgeo_RPANDA.sig2),sqrt(tf$DDexpgeo_mvMORPH.sig2))
abline(a=0,b=1,col="red")
plot(tf$DDexpgeo_RPANDA.z0,tf$DDexpgeo_mvMORPH.z0)
abline(a=0,b=1,col="red")

tf<-read.csv("~/ddexp/tests/20180802/tanager_fits_fruit_TESTS_20180804.csv")

par(mfrow=c(2,2))
plot(tf$DDexpgeo_RPANDA.lnL,tf$DDexpgeo_mvMORPH.lnL)
abline(a=0,b=1,col="red")
plot(tf$DDexpgeo_RPANDA.r,tf$DDexpgeo_mvMORPH.r)
abline(a=0,b=1,col="red")
plot(sqrt(tf$DDexpgeo_RPANDA.sig2),sqrt(tf$DDexpgeo_mvMORPH.sig2))
abline(a=0,b=1,col="red")
plot(tf$DDexpgeo_RPANDA.z0,tf$DDexpgeo_mvMORPH.z0)
abline(a=0,b=1,col="red")

max(abs(tf[,10]-tf[,16]))
#Copyright 2023 Tuobang Li

#These codes and manuscripts are under review in PNAS, please do not share them.

#If you are interested, please do not hesitate to contact me. Cooperation is also welcomed!

#require foreach and doparallel for parallel processing of bootstrap (not available for some types of computers)
if (!require("foreach")) install.packages("foreach")
library(foreach)
if (!require("doParallel")) install.packages("doParallel")
library(doParallel)
#require randtoolbox for random number generations
if (!require("randtoolbox")) install.packages("randtoolbox")
library(randtoolbox)
if (!require("Rcpp")) install.packages("Rcpp")
library(Rcpp)
if (!require("Rfast")) install.packages("Rfast")
library(Rfast)
if (!require("REDSReview")) install.packages("REDSReview_1.0.tar.gz", repos = NULL)
library(REDSReview)
if (!require("matrixStats")) install.packages("matrixStats")
library(matrixStats)
if (!require("rootSolve")) install.packages("rootSolve")
library(rootSolve)
numCores <- detectCores()-4 # Detect the number of available cores
cl <- makeCluster(numCores) # Create a cluster with the number of cores
registerDoParallel(cl) # Register the parallel backend


#set the stop criterion
criterionset=1e-10

#bootsize for bootstrap approximation of the distributions of the kernel of U-statistics.
n <- 13824*2*3
(n%%10)==0
# maximum order of moments
morder <- 6
#large sample size (approximating asymptotic)
largesize1<-13824*2

#generate quasirandom numbers based on the Sobol sequence
quasiunisobol<-sobol(n=n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                     mixed = FALSE, method = "C", start = 1)

quasiuni<-quasiunisobol

quasiuni_sorted2 <- na.omit(rowSort(quasiuni[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted3 <- na.omit(rowSort(quasiuni[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted4 <- na.omit(rowSort(quasiuni[,1:4], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted5 <- na.omit(rowSort(quasiuni[,1:5], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted6 <- na.omit(rowSort(quasiuni[,1:6], descend = FALSE, stable = FALSE, parallel = TRUE))


# Forever...

largesize1<-13824*2
samplesize=576*9
batchsizebase=1000
orderlist1_AB20<-createorderlist(quni1=quasiuni_sorted2,size=samplesize,interval=8,dimension=2)
orderlist1_AB20<-orderlist1_AB20[1:largesize1,]
orderlist1_AB30<-createorderlist(quni1=quasiuni_sorted3,size=samplesize,interval=8,dimension=3)
orderlist1_AB30<-orderlist1_AB30[1:largesize1,]
orderlist1_AB40<-createorderlist(quni1=quasiuni_sorted4,size=samplesize,interval=8,dimension=4)
orderlist1_AB40<-orderlist1_AB40[1:largesize1,]
orderlist1_AB50<-createorderlist(quni1=quasiuni_sorted5,size=samplesize,interval=8,dimension=5)
orderlist1_AB50<-orderlist1_AB50[1:largesize1,]
orderlist1_AB60<-createorderlist(quni1=quasiuni_sorted6,size=samplesize,interval=8,dimension=6)
orderlist1_AB60<-orderlist1_AB60[1:largesize1,]

morder=6

quasiunisobol0<-matrix(randtoolbox::SFMT(largesize1*3*morder),ncol=morder)

quasiuni0<-rbind(quasiunisobol0)

quasiunisobol0<-c()

quasiuni_sorted20 <- na.omit(rowSort(quasiuni0[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted30 <- na.omit(rowSort(quasiuni0[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted40 <- na.omit(rowSort(quasiuni0[,1:4], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted50 <- na.omit(rowSort(quasiuni0[,1:5], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted60 <- na.omit(rowSort(quasiuni0[,1:6], descend = FALSE, stable = FALSE, parallel = TRUE))

orderlist1_AB20_rand <-createorderlist(quni1=quasiuni_sorted20 ,size=samplesize,interval=16,dimension=2)
orderlist1_AB20_rand <-orderlist1_AB20_rand [1:largesize1,]
orderlist1_AB3_rand <-createorderlist(quni1=quasiuni_sorted30 ,size=samplesize,interval=16,dimension=3)
orderlist1_AB3_rand <-orderlist1_AB3_rand [1:largesize1,]
orderlist1_AB4_rand <-createorderlist(quni1=quasiuni_sorted40 ,size=samplesize,interval=16,dimension=4)
orderlist1_AB4_rand <-orderlist1_AB4_rand [1:largesize1,]
orderlist1_AB5_rand <-createorderlist(quni1=quasiuni_sorted50 ,size=samplesize,interval=16,dimension=5)
orderlist1_AB5_rand <-orderlist1_AB5_rand [1:largesize1,]
orderlist1_AB6_rand <-createorderlist(quni1=quasiuni_sorted60 ,size=samplesize,interval=16,dimension=6)
orderlist1_AB6_rand <-orderlist1_AB6_rand [1:largesize1,]
quasiuni_sorted20 <-c()
quasiuni_sorted30 <-c()
quasiuni_sorted40 <-c()
quasiuni_sorted50 <-c()
quasiuni_sorted60 <-c()



setSeed(1)
morder=6
quasiuni_M<-sobol(n=(largesize1*3*morder), dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                  mixed = FALSE, method = "C", start = 1)
largesize1<-13824*2
samplesize=576*9
orderlist1_hlsmall<-createorderlist(quni1=quasiuni_M[,1:6],size=samplesize,interval=8,dimension=6)
orderlist1_hlsmall<-orderlist1_hlsmall[1:largesize1,]
orderlist1_hllarge<-createorderlist(quni1=quasiuni_M[,1:6],size=largesize1,interval=8,dimension=6)
orderlist1_hllarge<-orderlist1_hllarge[1:largesize1,]

morder=6
# unibatchran_M<-matrix(randtoolbox::SFMT(largesize1*3*morder),ncol=morder)
# 
# orderlist1_hlsmall_rand<-createorderlist(quni1=unibatchran_M[,1:6],size=samplesize,interval=8,dimension=6)
# orderlist1_hlsmall_rand<-orderlist1_hlsmall_rand[1:largesize1,]

setSeed(1)

orderlist1_AB6_randomall<-c()
for(i in (1:batchsizebase)){
  unibatchran_M<-matrix(randtoolbox::SFMT(largesize1*3*morder),ncol=morder)
  
  orderlist1_AB6_random<-createorderlist(quni1=unibatchran_M[,1:6],size=samplesize,interval=8,dimension=6)
  orderlist1_AB6_random<-orderlist1_AB6_random[1:largesize1,]
  
  unibatchran_M<-c()
  
  orderlist1_AB6_randomall<-cbind(orderlist1_AB6_randomall,orderlist1_AB6_random)
}


batchsize=batchsizebase

n <- samplesize
setSeed(1)
unibatchran<-matrix(SFMT(samplesize*batchsize),ncol=batchsize)

unibatch<-colSort(unibatchran, descend = FALSE, stable = FALSE, parallel = TRUE)


kurtPareto<- read.csv(("kurtPareto_91210.csv"))
allkurtPareto<-unlist(kurtPareto)



simulatedbatch_ABSE<-foreach(batchnumber =c((1:length(allkurtPareto))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(randtoolbox)
  library(REDSReview)
  setSeed(1)
  set.seed(1)
  
  a=allkurtPareto[batchnumber]
  targetm<-a/(a-1)
  targetvar<-(((a))*(1)/((-2+(a))*((-1+(a))^2)))
  targettm<-((((a)+1)*(2)*(sqrt(a-2)))/((-3+(a))*(((a))^(1/2))))*(((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^3))
  targetfm<-(3+(6*((a)^3+(a)^2-6*(a)-2)/(((a))*((-3+(a)))*((-4+(a))))))*((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  
  SEbataches<-c()
  for (batch1 in c(1:batchsize)){
    x<-c(dsPareto(uni=unibatch[,batch1], shape=a, scale = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    
    Huberx<-Huber_estimator(x=sortedx, tol = 1e-10)
    SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted=TRUE)
    
    midhinge1<-midhinge(x=sortedx,sorted = TRUE)
    SWA81<-SWA8(x=sortedx,interval=8,batch="auto",sorted = TRUE)
    SWAHlmean1<-SWAHLmean(x=sortedx,orderlist1_sorted2=orderlist1_AB20,orderlist1_sorted3=orderlist1_AB30,orderlist1_sorted4=orderlist1_AB40,orderlist1_sorted5=orderlist1_AB50,orderlist1_sorted6=orderlist1_AB60,batch="auto")
    
    SWAHlmean1_rand<-SWAHLmean(x=sortedx,orderlist1_sorted2=orderlist1_AB20_rand,orderlist1_sorted3=orderlist1_AB3_rand,orderlist1_sorted4=orderlist1_AB4_rand ,orderlist1_sorted5=orderlist1_AB5_rand ,orderlist1_sorted6=orderlist1_AB6_rand ,batch="auto")
    
    rqm1<-rqm(x=sortedx,sorted = TRUE)
    
    MoM2<-median_of_means(sortedx,korder=2)
    MoM3<-median_of_means(sortedx,korder=3)
    MoM4<-median_of_means(sortedx,korder=4)
    MoM5<-median_of_means(sortedx,korder=5)
    MoM519<-median_of_means(sortedx,korder=log(0.5,7/8))
    
    Hodges_Lehmann_estimator<-wilcox.test(x=sortedx,conf.int = TRUE)$estimate
    SWA22<-SWA(x = sortedx, percentage = 1 - (1 -1/2)^0.5, batch = "auto", sorted = TRUE, rand = TRUE)
    
    
    sortedx<-c()

    allrawmoBias<-c(
      firstbias=(c(Huberx,Hodges_Lehmann_estimator=Hodges_Lehmann_estimator,SWA22,SMWM9,midhinge1,SWA81=SWA81,SWAHlmean1=SWAHlmean1,SWAHlmean1_rand=SWAHlmean1_rand,rqm1,MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM519=MoM519)-targetm)
    )
    
    allrawmo1<-c(first=(c(Huberx,Hodges_Lehmann_estimator=Hodges_Lehmann_estimator,SWA22,SMWM9,midhinge1,SWA81=SWA81,SWAHlmean1=SWAHlmean1,SWAHlmean1_rand=SWAHlmean1_rand,rqm1,MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoM519=MoM519))
    )
    
    all1<-(c(kurtx=kurtx,skewx=skewx,allrawmoBias,allrawmo1,targetall))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  
  write.csv(SEbataches,paste("Pareto_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,3:82])^2))/sqrt(targetvar)
  
  
  AB1_mean<-abs(colMeans((SEbataches[,3:82])))/sqrt(targetvar)
  
  
  SEbatachesmean <- colMeans(SEbataches)
  
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,83:166]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  
  allSE<-c(mean_SE1=mean_SE1)
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1  )

  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],RMSE1_mean=RMSE1_mean,AB1_mean=AB1_mean,allSE=allSE,allSE_unstan=allSE_unstan,SEbatachesmean=SEbatachesmean)
  
  
  allErrors
}


write.csv(simulatedbatch_ABSE,paste("Pareto_ABSSE.csv", sep = ","), row.names = FALSE)


simulatedbatch_ABSE_SE<-foreach(batchnumber =c((1:length(allkurtPareto))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(randtoolbox)
  library(REDSReview)
  setSeed(1)
  set.seed(1)
  
  
  a=allkurtPareto[batchnumber]
  targetm<-a/(a-1)
  targetvar<-(((a))*(1)/((-2+(a))*((-1+(a))^2)))
  targettm<-((((a)+1)*(2)*(sqrt(a-2)))/((-3+(a))*(((a))^(1/2))))*(((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^3))
  targetfm<-(3+(6*((a)^3+(a)^2-6*(a)-2)/(((a))*((-3+(a)))*((-4+(a))))))*((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  
  SEbataches<- read.csv(paste("Pareto_raw_ABSSE_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  
  SEbatachesmean <- colMeans(SEbataches)
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,83:166]), 2, se_mean)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)

  
  allSE<-c(mean_SE1=mean_SE1,SEbatachesmean[1]  )
  allSE_unstan<-c(SEbatachesmean[1],meansd_unscaled1=meansd_unscaled1  )
  
  
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],skew=SEbatachesmean[2],se_mean=se_mean_all1,allSE=allSE,allSE_unstan=allSE_unstan,SEbatachesmean=SEbatachesmean)
  
  allErrors
}



write.csv(simulatedbatch_ABSE_SE,paste("Pareto_ABSSE_error.csv", sep = ","), row.names = FALSE)




stopCluster(cl)
registerDoSEQ()


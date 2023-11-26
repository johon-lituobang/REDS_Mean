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


numCores <- detectCores()-4
#registering clusters, can set a smaller number using numCores-1

registerDoParallel(numCores)


#bootsize for bootstrap approximation of the distributions of the kernel of U-statistics.
n <- 2048*9*3
(n%%10)==0
# maximum order of moments
morder <- 6

#generate quasirandom numbers based on the Sobol sequence
quasiunisobol<-sobol(n=n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                     mixed = FALSE, method = "C", start = 1)
quasiuni<-rbind(quasiunisobol)

quasiunisobol<-c()

quasiuni_sorted2 <- na.omit(rowSort(quasiuni[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted3 <- na.omit(rowSort(quasiuni[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted4 <- na.omit(rowSort(quasiuni[,1:4], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted5 <- na.omit(rowSort(quasiuni[,1:5], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted6 <- na.omit(rowSort(quasiuni, descend = FALSE, stable = FALSE, parallel = TRUE))

# Forever...

#set the stop criterion
criterionset=1e-30




kurtgnorm<- read.csv(("kurtgnorm_31150.csv"))
allkurtgnorm<-unlist(kurtgnorm)

largesize1=2048*9
samplesize=2048*2
batchsizebase=1000
orderlist1_AB2<-createorderlist(quni1=quasiuni_sorted2,size=samplesize,interval=8,dimension=2)
orderlist1_AB2<-orderlist1_AB2[1:largesize1,]
orderlist1_AB3<-createorderlist(quni1=quasiuni_sorted3,size=samplesize,interval=8,dimension=3)
orderlist1_AB3<-orderlist1_AB3[1:largesize1,]
orderlist1_AB4<-createorderlist(quni1=quasiuni_sorted4,size=samplesize,interval=8,dimension=4)
orderlist1_AB4<-orderlist1_AB4[1:largesize1,]
orderlist1_AB5<-createorderlist(quni1=quasiuni_sorted5,size=samplesize,interval=8,dimension=5)
orderlist1_AB5<-orderlist1_AB5[1:largesize1,]
orderlist1_AB6<-createorderlist(quni1=quasiuni_sorted6,size=samplesize,interval=8,dimension=6)
orderlist1_AB6<-orderlist1_AB6[1:largesize1,]

batchsize=batchsizebase

n <- samplesize
setSeed(1)
unibatchran<-matrix(SFMT(samplesize*batchsize),ncol=batchsize)

unibatch<-colSort(unibatchran, descend = FALSE, stable = FALSE, parallel = TRUE)


simulatedbatch_ABSE<-foreach(batchnumber =c((1:length(allkurtgnorm))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(randtoolbox)
  library(REDSReview)
  library(rootSolve)
  setSeed(1)
  set.seed(1)

  a=allkurtgnorm[batchnumber]
  
  targetm<-0
  targetvar<-gamma(3/a)/((gamma(1/a)))
  targettm<-0
  targetfm<-((gamma(3/a)/((gamma(1/a))))^2)*gamma(5/a)*gamma(1/a)/((gamma(3/a))^2)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  SEbataches<-c()
  for (batch1 in c(1:batchsize)){
    x<-c(dsgnorm(uni=unibatch[,batch1], shape=a/1, scale = 1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
    x<-c()
    Huberx<-Huber_estimator(x=sortedx)
    SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted = TRUE)
    midhinge1<-midhinge(x=sortedx,sorted = TRUE)
    SWA81<-SWA8(x=sortedx,interval=8,batch="auto",sorted = TRUE)
    SWAHlmean1<-SWAHLmean(x=sortedx,orderlist1_sorted2=orderlist1_AB2,orderlist1_sorted3=orderlist1_AB3,orderlist1_sorted4=orderlist1_AB4,orderlist1_sorted5=orderlist1_AB5,orderlist1_sorted6=orderlist1_AB6,batch="auto")
    SWA22<-SWA(x = sortedx, percentage = 1 - (1 -1/2)^0.5, batch = "auto", sorted = TRUE, rand = TRUE)
    
    rqm1<-rqm(x=sortedx,sorted = TRUE)
    MoM2<-median_of_means(sortedx,korder=2)
    MoM3<-median_of_means(sortedx,korder=3)
    MoM4<-median_of_means(sortedx,korder=4)
    MoM5<-median_of_means(sortedx,korder=5)
    MoRM2<-mHLM(x=sortedx,dimension=2,boot=TRUE,quasi=FALSE,largesize=largesize1,interval=8)
    MoRM3<-mHLM(x=sortedx,dimension=3,boot=TRUE,quasi=FALSE,largesize=largesize1,interval=8)
    MoRM4<-mHLM(x=sortedx,dimension=4,boot=TRUE,quasi=FALSE,largesize=largesize1,interval=8)
    MoRM5<-mHLM(x=sortedx,dimension=5,boot=TRUE,quasi=FALSE,largesize=largesize1,interval=8)
    momentsx<-unbiasedmoments(x=sortedx)
    sortedx<-c()
    
    allrawmoBias<-c(
      firstbias=(c(Huberx,SMWM9,midhinge1,SWA81,SWAHlmean1,SWAHL=SWA22,rqm1,MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoRM2=MoRM2,MoRM3=MoRM3,MoRM4=MoRM4,MoRM5=MoRM5)-targetm)
    )
    
    all1<-t(c(kurtx,Huberx,SMWM9,midhinge1,SWA81,SWAHlmean1,SWAHL=SWA22,rqm1,MoM2=MoM2,MoM3=MoM3,MoM4=MoM4,MoM5=MoM5,MoRM2=MoRM2,MoRM3=MoRM3,MoRM4=MoRM4,MoRM5=MoRM5,targetall,momentsx,allrawmoBias))
    
    SEbataches<-rbind(SEbataches,all1)
  }
  
  write.csv(SEbataches,paste("gnorm_raw_ABSE_SRM_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","), row.names = FALSE)
  
  RMSE1_mean<-sqrt(colMeans((SEbataches[,64:117])^2))/sqrt(targetvar)
  
  AB1_mean<-abs(colMeans((SEbataches[,64:117])))/sqrt(targetvar)
  
  SEbatachesmean <- colMeans(SEbataches)
  
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,c(1:63)]), 2, unbiasedsd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  # mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,c(1:60)])/ratiomean1)
  # 
  # meansd1<-apply((mean_SEbatachesmeanprocess), 2, unbiasedsd)
  # mean_SSE1<-meansd1/sqrt(targetvar)
  # 
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],RMSE1_mean=RMSE1_mean,AB1_mean=AB1_mean,SEbatachesmean=SEbatachesmean,meansd_unscaled1=meansd_unscaled1,mean_SE1=mean_SE1)

  allErrors
}

write.csv(simulatedbatch_ABSE,paste("gnorm_ABSE_SRM.csv", sep = ","), row.names = FALSE)



simulatedbatch_ABSE_SE<-foreach(batchnumber =c((1:length(allkurtgnorm))), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(randtoolbox)
  library(REDSReview)
  library(rootSolve)
  setSeed(1)
  set.seed(1)
  
  
  a=allkurtgnorm[batchnumber]
  
  targetm<-0
  targetvar<-gamma(3/a)/((gamma(1/a)))
  targettm<-0
  targetfm<-((gamma(3/a)/((gamma(1/a))))^2)*gamma(5/a)*gamma(1/a)/((gamma(3/a))^2)
  kurtx<-targetfm/(targetvar^(4/2))
  skewx<-targettm/(targetvar^(3/2))
  
  SEbataches<- read.csv(paste("gnorm_raw_ABSE_SRM_finite",samplesize,round(kurtx,digits = 1),".csv", sep = ","))
  
  
  SEbatachesmean <- colMeans(SEbataches)
  
  
  meansd_unscaled1<-apply((SEbataches[1:batchsize,c(1:63)]), 2, se_sd)
  
  mean_SE1<-meansd_unscaled1/sqrt(targetvar)
  # mean_SEbatachesmeanprocess<-t(t(SEbataches[1:batchsize,c(1:60)])/ratiomean1)
  # 
  # meansd1<-apply((mean_SEbatachesmeanprocess), 2, se_sd)
  # mean_SSE1<-meansd1/sqrt(targetvar)
  # 
  se_mean_all1<-apply((SEbataches[1:batchsize,]), 2, se_mean)
  
  allErrors<-c(samplesize=samplesize,kurt=SEbatachesmean[1],se_mean=se_mean_all1,SEbatachesmean=SEbatachesmean,meansd_unscaled1=meansd_unscaled1,mean_SE1=mean_SE1)

  allErrors
}

write.csv(simulatedbatch_ABSE_SE,paste("gnorm_ABSE_SRM_error.csv", sep = ","), row.names = FALSE)


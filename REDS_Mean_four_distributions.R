#Copyright 2023 Tuobang Li

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
numCores <- detectCores()-4
#registering clusters, can set a smaller number using numCores-1

registerDoParallel(numCores)


asymptotic_n <- 2048*900*2*3
(asymptotic_n%%10)==0
# maximum order of kernel
morder <- 6
#large sample size (asymptotic bias)
largesize<-2048*900*2

#generate quasirandom numbers based on the Sobol sequence
quasiunisobol_asymptotic<-sobol(n=asymptotic_n, dim = morder, init = TRUE, scrambling = 0, seed = NULL, normal = FALSE,
                                mixed = FALSE, method = "C", start = 1)

quasiuni_asymptotic<-rbind(quasiunisobol_asymptotic)

quasiunisobol_asymptotic<-c()

quasiuni_sorted2_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:2], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted3_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:3], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted4_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:4], descend = FALSE, stable = FALSE, parallel = TRUE))

quasiuni_sorted5_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:5], descend = FALSE, stable = FALSE, parallel = TRUE))
quasiuni_sorted6_asymptotic <- na.omit(rowSort(quasiuni_asymptotic[,1:6], descend = FALSE, stable = FALSE, parallel = TRUE))

# Forever...

quasiuni_asymptotic<-(quasiuni_asymptotic[,1])
quasiuni_asymptotic<-quasiuni_asymptotic[1:largesize]
orderlist1_AB2_asymptotic<-createorderlist(quni1=quasiuni_sorted2_asymptotic,size=largesize,interval=8,dimension=2)
orderlist1_AB2_asymptotic<-orderlist1_AB2_asymptotic[1:largesize,]
orderlist1_AB3_asymptotic<-createorderlist(quni1=quasiuni_sorted3_asymptotic,size=largesize,interval=8,dimension=3)
orderlist1_AB3_asymptotic<-orderlist1_AB3_asymptotic[1:largesize,]
orderlist1_AB4_asymptotic<-createorderlist(quni1=quasiuni_sorted4_asymptotic,size=largesize,interval=8,dimension=4)
orderlist1_AB4_asymptotic<-orderlist1_AB4_asymptotic[1:largesize,]
orderlist1_AB5_asymptotic<-createorderlist(quni1=quasiuni_sorted5_asymptotic,size=largesize,interval=8,dimension=5)
orderlist1_AB5_asymptotic<-orderlist1_AB5_asymptotic[1:largesize,]
orderlist1_AB6_asymptotic<-createorderlist(quni1=quasiuni_sorted6_asymptotic,size=largesize,interval=8,dimension=6)
orderlist1_AB6_asymptotic<-orderlist1_AB6_asymptotic[1:largesize,]

quasiuni_sorted2_asymptotic<-c()
quasiuni_sorted3_asymptotic<-c()
quasiuni_sorted4_asymptotic<-c()
quasiuni_sorted5_asymptotic<-c()
quasiuni_sorted6_asymptotic<-c()
#These codes and manuscripts are under review in PNAS, please do not share them.

#If you are interested, please do not hesitate to contact me. Coorperation is also welcomed.

#Weibull

kurtWeibull<- read.csv(("kurtWeibull_31150.csv"))
allkurtWeibull<-unlist(kurtWeibull)

simulatedbatchWeibull_asymptoticbias<-foreach(batchnumber = (1:length(allkurtWeibull)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(REDSReview)
  library(rootSolve)
  a=allkurtWeibull[batchnumber]
  x<-c(dsWeibull(uni=quasiuni_asymptotic, shape=a/1, scale = 1))
  targetm<-gamma(1+1/(a/1))
  targetvar<-(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2)
  targettm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^3)*(gamma(1+3/(a/1))-3*(gamma(1+1/(a/1)))*((gamma(1+2/(a/1))))+2*((gamma(1+1/(a/1)))^3))/((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(3))
  targetfm<-((sqrt(gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^4)*(gamma(1+4/(a/1))-4*(gamma(1+3/(a/1)))*((gamma(1+1/(a/1))))+6*(gamma(1+2/(a/1)))*((gamma(1+1/(a/1)))^2)-3*((gamma(1+1/(a/1)))^4))/(((gamma(1+2/(a/1))-(gamma(((1+1/(a/1)))))^2))^(2))
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  
  targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
  x<-c()
  Huberx<-Huber_estimator(x=sortedx)
  SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted = TRUE)
  midhinge1<-midhinge(x=sortedx,sorted = TRUE)
  SWA81<-SWA8(x=sortedx,interval=8,batch="auto",sorted = TRUE)
  SWAHlmean1<-SWAHLmean(x=sortedx,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,orderlist1_sorted5=orderlist1_AB5_asymptotic,orderlist1_sorted6=orderlist1_AB6_asymptotic,batch="auto")
  SWA22<-SWA(x = sortedx, percentage =1 - (1 - 1/2)^0.5, batch = "auto", sorted = TRUE, rand = TRUE)
  rqm1<-rqm(x=sortedx,sorted = TRUE)
  momentsx<-unbiasedmoments(x=sortedx)
  sortedx<-c()
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,midhinge1,SWA81,SWAHlmean1,SWAHL=SWA22,rqm1)-targetm)/(sqrt(targetvar))
  )
  
  all1<-(c(kurtx,Huberx,SMWM9,midhinge1,SWA81,SWAHlmean1,SWAHL=SWA22,rqm1,targetall,momentsx,allrawmoBias))
}

write.csv(simulatedbatchWeibull_asymptoticbias,paste("asymptotic_Weibull_SRM_Process",largesize,".csv", sep = ","), row.names = FALSE)

#gamma

kurtgamma<- read.csv(("kurtgamma_31150.csv"))
allkurtgamma<-unlist(kurtgamma)

simulatedbatchgamma_asymptoticbias<-foreach(batchnumber = (1:length(allkurtgamma)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(REDSReview)
  library(rootSolve)
  a=allkurtgamma[batchnumber]
  x<-c(dsgamma(uni=quasiuni_asymptotic, shape=a, rate   = 1))
  targetm<-a
  targetvar<-(a)
  targettm<-((sqrt(a))^3)*2/sqrt(a)
  targetfm<-((sqrt(a))^4)*((6/(a))+3)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  
  targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
  x<-c()
  Huberx<-Huber_estimator(x=sortedx)
  SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted = TRUE)
  midhinge1<-midhinge(x=sortedx,sorted = TRUE)
  SWA81<-SWA8(x=sortedx,interval=8,batch="auto",sorted = TRUE)
  SWAHlmean1<-SWAHLmean(x=sortedx,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,orderlist1_sorted5=orderlist1_AB5_asymptotic,orderlist1_sorted6=orderlist1_AB6_asymptotic,batch="auto")
  SWA22<-SWA(x = sortedx, percentage =1 - (1 - 1/2)^0.5, batch = "auto", sorted = TRUE, rand = TRUE)
  rqm1<-rqm(x=sortedx,sorted = TRUE)
  momentsx<-unbiasedmoments(x=sortedx)
  sortedx<-c()
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,midhinge1,SWA81,SWAHlmean1,SWAHL=SWA22,rqm1)-targetm)/(sqrt(targetvar))
  )
  
  all1<-(c(kurtx,Huberx,SMWM9,midhinge1,SWA81,SWAHlmean1,SWAHL=SWA22,rqm1,targetall,momentsx,allrawmoBias))
}

write.csv(simulatedbatchgamma_asymptoticbias,paste("asymptotic_gamma_SRM_Process",largesize,".csv", sep = ","), row.names = FALSE)


#Pareto
kurtPareto<- read.csv(("kurtPareto_91210.csv"))
allkurtPareto<-unlist(kurtPareto)

simulatedbatchPareto_asymptoticbias<-foreach(batchnumber = (1:length(allkurtPareto)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(REDSReview)
  library(rootSolve)
  a=allkurtPareto[batchnumber]
  x<-c(dsPareto(uni=quasiuni_asymptotic, shape=a, scale = 1))
  targetm<-a/(a-1)
  targetvar<-(((a))*(1)/((-2+(a))*((-1+(a))^2)))
  targettm<-((((a)+1)*(2)*(sqrt(a-2)))/((-3+(a))*(((a))^(1/2))))*(((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^3))
  targetfm<-(3+(6*((a)^3+(a)^2-6*(a)-2)/(((a))*((-3+(a)))*((-4+(a))))))*((sqrt(((a))*(1)/((-2+(a))*((-1+(a))^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  
  targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
  x<-c()
  Huberx<-Huber_estimator(x=sortedx)
  SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted = TRUE)
  midhinge1<-midhinge(x=sortedx,sorted = TRUE)
  SWA81<-SWA8(x=sortedx,interval=8,batch="auto",sorted = TRUE)
  SWAHlmean1<-SWAHLmean(x=sortedx,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,orderlist1_sorted5=orderlist1_AB5_asymptotic,orderlist1_sorted6=orderlist1_AB6_asymptotic,batch="auto")
  SWA22<-SWA(x = sortedx, percentage =1 - (1 - 1/2)^0.5, batch = "auto", sorted = TRUE, rand = TRUE)
  rqm1<-rqm(x=sortedx,sorted = TRUE)
  momentsx<-unbiasedmoments(x=sortedx)
  sortedx<-c()
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,midhinge1,SWA81,SWAHlmean1,SWAHL=SWA22,rqm1)-targetm)/(sqrt(targetvar))
  )
  
  all1<-(c(kurtx,Huberx,SMWM9,midhinge1,SWA81,SWAHlmean1,SWAHL=SWA22,rqm1,targetall,momentsx,allrawmoBias))
}

write.csv(simulatedbatchPareto_asymptoticbias,paste("asymptotic_Pareto_SRM_Process",largesize,".csv", sep = ","), row.names = FALSE)



#lognorm

kurtlognorm<- read.csv(("kurtlognorm_31150.csv"))
allkurtlognorm<-unlist(kurtlognorm)

simulatedbatchlognorm_asymptoticbias<-foreach(batchnumber = (1:length(allkurtlognorm)), .combine = 'rbind') %dopar% {
  library(Rfast)
  library(matrixStats)
  library(REDSReview)
  library(rootSolve)
  a=allkurtlognorm[batchnumber]
  x<-c(dslnorm(uni=quasiuni_asymptotic, meanlog =0, sdlog   = a/1))
  targetm<-exp((a^2)/2)
  targetvar<-(exp((a/1)^2)*(-1+exp((a/1)^2)))
  targettm<-sqrt(exp((a/1)^2)-1)*((2+exp((a/1)^2)))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^3)
  targetfm<-((-3+exp(4*((a/1)^2))+2*exp(3*((a/1)^2))+3*exp(2*((a/1)^2))))*((sqrt(exp((a/1)^2)*(-1+exp((a/1)^2))))^4)
  kurtx<-targetfm/(targetvar^(4/2))
  kurtx<-c(kurtx=kurtx)
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  
  targetall<-c(targetm=targetm,targetvar=targetvar,targettm=targettm,targetfm=targetfm)
  x<-c()
  Huberx<-Huber_estimator(x=sortedx)
  SMWM9<-SWA9(x=sortedx,interval=9,batch="auto",sorted = TRUE)
  midhinge1<-midhinge(x=sortedx,sorted = TRUE)
  SWA81<-SWA8(x=sortedx,interval=8,batch="auto",sorted = TRUE)
  SWAHlmean1<-SWAHLmean(x=sortedx,orderlist1_sorted2=orderlist1_AB2_asymptotic,orderlist1_sorted3=orderlist1_AB3_asymptotic,orderlist1_sorted4=orderlist1_AB4_asymptotic,orderlist1_sorted5=orderlist1_AB5_asymptotic,orderlist1_sorted6=orderlist1_AB6_asymptotic,batch="auto")
  SWA22<-SWA(x = sortedx, percentage =1 - (1 - 1/2)^0.5, batch = "auto", sorted = TRUE, rand = TRUE)
  rqm1<-rqm(x=sortedx,sorted = TRUE)
  momentsx<-unbiasedmoments(x=sortedx)
  sortedx<-c()
  allrawmoBias<-c(
    firstbias=abs(c(Huberx,SMWM9,midhinge1,SWA81,SWAHlmean1,SWAHL=SWA22,rqm1)-targetm)/(sqrt(targetvar))
  )
  
  all1<-(c(kurtx,Huberx,SMWM9,midhinge1,SWA81,SWAHlmean1,SWAHL=SWA22,rqm1,targetall,momentsx,allrawmoBias))
}

write.csv(simulatedbatchlognorm_asymptoticbias,paste("asymptotic_lognorm_SRM_Process",largesize,".csv", sep = ","), row.names = FALSE)



registerDoSEQ()


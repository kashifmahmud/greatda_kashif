library("plantecophys")
library("doBy")
library("nlme")
library("dplyr")
library(readxl)
library(stringr)
library("nlstools")
library(plotBy)
library(magicaxis)
library(rgl)
library( propagate)
library(gplots)
library(ggplot2)
library(lme4)
library(lubridate)
library(car)
library(boot)
library(multcomp)
library(plyr)
library(knitr)
library(chron)
library(sp)
library(devtools)
library(HIEv)
library(dismo)
library(rgbif)  # when using GBIF
library(ALA4R)  # when using ALA
library(raster)
library(oz)
library(installr)
library(broom)
library(piecewiseSEM)
library(plover)
#-----------------------------------------------------------------------------------------------------------------------
#function to read and process wtcIII ACi data 

getwtc3<-function(path){
  #read wtcIII ACi data
  
  aci<-read.csv(paste(path,"/WTC_TEMP_CM_GX-ACITEMP_20130917-20140507_R_V2.csv",sep=""))
  
  
  
  #to define months
  aci$Season<-NA
  aci$Season[which(aci$MeanTemp=="14-Jan")]<-"January"
  aci$Season[which(aci$MeanTemp=="13-Sep")]<-"September"
  aci$Season[which(aci$MeanTemp=="14-Apr")]<-"April"
  aci$Season[which(aci$MeanTemp=="14-May")]<-"April"
  aci$Season<-as.factor(aci$Season)
  
  #to drop TargetTemp 35. 
  #stomatal conductance of this temperature level is very low (less than 0.09 for all 31 measurements)
  
  #aci<-aci[!aci$TargetTemp==35,]
  aci$Curve<-aci$Plotsite
  #write.csv(aci,"C:/Dushan/Repos//ACI_ALL/WTC3/WTC3_ACidata_processed.csv",row.names=FALSE,sep=",")
  
  return(aci)  
}


#------------------------------------------------------------------------------------------------------------------------

#function to read and process wtcII ACi data 

getwtc2<-function(path){
  #read wtcIII ACi data
  
  aciG<-read.csv(paste(path,"/WTC_Temp&CO2_CM_GX-AciTemp_20101202-20110908_R.csv",sep=""))
  
  #to define temperatures
  aciG$TargetTemp<-NA
  aciG$TargetTemp[which(aciG$Species=="Eglob @ 25C")]<-25
  aciG$TargetTemp[which(aciG$Species=="Eglob @ 20C")]<-20
  aciG$TargetTemp[which(aciG$Species=="Eglob @ 30C")]<-30
  aciG$TargetTemp[which(aciG$Species=="Eglob @ 35C")]<-35
  aciG$TargetTemp[which(aciG$Species=="Eglob @ 41C")]<-41
  aciG$TargetTemp[which(aciG$Species=="Eglob @ 17C")]<-17
  aciG$TargetTemp[which(aciG$Species=="Eglob")]<-25
  aciG$TargetTemp<-as.numeric(aciG$TargetTemp)
  
  #to extract month
  aciG$Date<- as.Date(mdy(aciG$Date))
  aciG$Month<-month(aciG$Date)
  
  
  
  
  aciG$Season<-ifelse(aciG$Month==02,"February",
                      ifelse(aciG$Month==08,"August",
                             ifelse(aciG$Month==09,"August",
                                    ifelse(aciG$Month==12,"December",
                                           ifelse(aciG$Month==04,"April",NA)))))
  
  aciG$Season<-as.factor(aciG$Season) 
  aciG$Curve<-aciG$Identity
  
  #write.table(aciG,"C:/Dushan/Repos/ACI_ALL/WTC2/WTC2_ACidata_processednew.csv",row.names=FALSE,sep=",")
  return(aciG)  
  

  
  
}


#------------------------------------------------------------------------------------------------------------------------
#function to read and process Tumbarumba ACi data 

getTb<-function(path="C:/Dushan/R/WTCIII"){
  
  #read  ACi data
  
  aciT<-read.csv(paste(path,"/rawdata/TumbarumbaGasex_ACis_Medlyn.csv",sep=""))
  
  #aciT$M<-month(as.POSIXlt(aciT$Date, format="%d/%m/%Y"))
  
  aciT$Date<- as.Date(dmy(aciT$Date))
  t<-data.frame(
    DateTime=aciT$Date,
    Season=format(as.POSIXct(aciT$Date, format="%d-%m-%y"), format="%m"))
  t$index<-c(1:702)
  aciT$index<-c(1:702)
  aciT<-merge(aciT,t,by=c("index"))
  
  aciT<-subset(aciT,aciT$PARi>1200)
  
  aciT$Species<-"E.delegatensis"

  aciT$Month<-NA
  aciT$Month[which(aciT$Season=="02")]<-"February"
  aciT$Month[which(aciT$Season=="05")]<-"May"
  aciT$Month[which(aciT$Season=="11")]<-"November"
  
  
  write.csv(aciT,"C:/Dushan/Repos//ACI_ALL/TMB/Tumbarumba_ACidata_processed.csv",row.names=FALSE)
  return(aciT)  
  
}
#getTb()
#-----------------------------------------------------------------------------------------------------------------------

#function to read and process WTCIII DIURNAL Photosynthesis data 
getD<-function(path="C:/Dushan/R/WTCIII"){
  
  #read  ACi data
  
  dataD<-read.csv(paste(path,"/rawdata/WTC_TEMP_CM_GX-DIURNAL_20130710-20140220_L1_v2.csv",sep=""))
  
  #Allocation of curve no for each canopy position in each campaign
  
  dataD$Number <- c(1)
  Number <- c()
  count <- 1  
  for (i in 2:length(dataD$campaign)){
    ifelse(dataD$position[i-1] == dataD$position[i],count <- count,count <- count + 1)
    Number[i] <- count 
  }
  dataD$Number[2:length(dataD$position)] <- na.omit(Number)
  
  return(dataD)
  
}

#-----------------------------------------------------------------------------------------------------------------------
#function to read and process HFE Commongarden ACi data (Lin's data) 

gethfe<-function(path){
  
  #read HFE Commongarden ACi data 
  
  acihfe<-read.csv(paste(path,"/HFE_CommonGarden_ACI_3Seasons_YSLin2012.csv",sep=""))
  

  #to define months
  acihfe$Month<-NA
  acihfe$Month[which(acihfe$season==1)]<-"August"
  acihfe$Month[which(acihfe$season==2)]<-"November"
  acihfe$Month[which(acihfe$season==3)]<-"February"
  
  acihfe$season<-as.factor(acihfe$season)
  
  acihfe$Spp<-NA
  acihfe$Spp[which(acihfe$Species=="Ecla")]<-"E.cladocalyx"
  acihfe$Spp[which(acihfe$Species=="Ecre")]<-"E.crebra"
  acihfe$Spp[which(acihfe$Species=="Edun")]<-"E.dunnii"
  acihfe$Spp[which(acihfe$Species=="Emel")]<-"E.melliodora"
  acihfe$Spp[which(acihfe$Species=="Esal")]<-"E.saligna"
  acihfe$Spp[which(acihfe$Species=="Emel")]<-"E.tereticornis"
  acihfe$Spp<-as.factor(acihfe$Spp)
  
  return(acihfe) 
  
  write.csv(acihfe,"C:/Dushan/Repos/PhotoM/Data/hfe_ACidata_processed.csv",row.names=FALSE)
}


#-----------------------------------------------------------------------------------------------------------------------


# fit temperature response of Vcmax and get CI of parameters


fVc <- as.formula(Vcmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)) * 
                    (1+exp((298.15*delS - 200)/(298.15*0.008314))) / 
                    (1+exp((TsK*delS-200)/(TsK*0.008314))))

fVJ <- as.formula(Jmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)) * 
                    (1+exp((298.15*delS - 200)/(298.15*0.008314))) / 
                    (1+exp((TsK*delS-200)/(TsK*0.008314))))


fVc.a<-as.formula(Vcmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)))
fVJ.a<-as.formula(Jmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)))

# peaked model (estimation of predictions and 95% CI dissabled to reduce impliment time)

fitpeaked<-function(dat){
  
  try(Vc<-nls(fVc, start = list(k25=100, Ea=60, delS = 0.64), data = dat))
  Vc2<-summary(Vc)
  res<-Vc2$coefficients[1:6]
  names(res)[1:6]<-c("Vcmax25","EaV","delsV","Vcmax25.se","EaV.se","delsV.se")
  
  #tt<-seq(min(dat$TsK),max(dat$TsK),length=51)
  #predicts<-predictNLS(Vc,newdata=data.frame(TsK=tt),interval="confidence",level=0.95)
 # predicts.df<-data.frame(predicts$summary)
  #predicts.df$TsK<-tt
  
  
  #return(list(results,predicts.df))
  r<-cor(fitted(Vc),dat$Vcmax)
  r2<-r*r
  
  #test for normality of residuals
  rest<-residuals(Vc)
  norm<-shapiro.test(rest)
  s<-norm$statistic
  pvalue<-norm$p.value
  
  topt<-200/(res[[3]]-0.008314*log(res[[2]]/(200-res[[2]])))
  Topt<-topt-273.15
  
  param<-c(res,Topt,r2,s,pvalue)
  names(param)[1:10]<-c("Vcmax25","Ea","delS","Vcmax25.se","Ea.se","delS.se","Topt","R2","S","pvalue")
  return(param)
  
#return(c(res,r2,s,pvalue))
  
  
  #return(list(res,predicts.df))
  #return(list(res))
}

# fit temperature response of Jmax and get CI of parameters
# peaked model (estimation of predictions and 95% CI dissabled to reduce impliment time)

fitpeakedJ<-function(dat){
  
  try(Vj<-nls(fVJ, start = list(k25=60, Ea=60, delS = 0.63), data = dat))
  Vj2<-summary(Vj)
  res<-Vj2$coefficients[1:6]
  names(res)[1:6]<-c("Jmax25","Ea","delsJ","Jmax25.se","EaJ.se","delsJ.se")
  
  #tt<-seq(min(dat$TsK),max(dat$TsK),length=51)
  #predicts<-predictNLS(Vj,newdata=data.frame(TsK=tt),interval="confidence",level=0.95)
  #predicts.df<-data.frame(predicts$summary)
  #predicts.df$TsK<-tt
  
  #return(list(results,predicts.df))
  r<-cor(fitted(Vj),dat$Jmax)
  r2<-r*r
  
  #test for normality of residuals
  rest<-residuals(Vj)
  norm<-shapiro.test(rest)
  s<-norm$statistic
  pvalue<-norm$p.value
  
  topt<-200/(res[[3]]-0.008314*log(res[[2]]/(200-res[[2]])))
  Topt<-topt-273.15
  
  param<-c(res,Topt,r2,s,pvalue)
  names(param)[1:10]<-c("Jmax25","Ea","delS","Jmax25.se","Ea.se","delS.se","Topt","R2","S","pvalue")
  return(param)
  
  #return(c(res,r2,s,pvalue))
  #return(list(res))
}

# function to fit temperature response of Vcmax and get CI for parameters
# standard Arrhenius model

fitarh<-function(dat){
  
  try(Vr<-nls(fVc.a, start = list(k25=100, Ea=40), data = dat))
  Vr2<-summary(Vr)
  res<-Vr2$coefficients[1:4]
  names(res)[1:4]<-c("Vcmax25","EaV","Vcmax25.se","EaV.se")
  
  #tt<-seq(min(dat$TsK),max(dat$TsK),length=51)
  #predicts<-predictNLS(Vr,newdata=data.frame(TsK=tt),interval="confidence",level=0.95)
  #predicts.df<-data.frame(predicts$summary)
  #predicts.df$TsK<-tt
  #return(list(results,predicts.df))
  
  r<-cor(fitted(Vr),dat$Vcmax)
  R2<-r*r
  
  #test for normality of residuals
  rest<-residuals(Vr)
  norm<-shapiro.test(rest)
  S<-norm$statistic
  pvalue<-norm$p.value
  
  param<-c(res,R2)
  names(param)[1:5]<-c("Vcmax25","EaV","Vcmax25.se","EaV.se","R2_Arr")
  return(param)
  
  
  #return(c(res,R2,S,pvalue))
  #return(list(res,predicts.df))
  #return(res)
}


#for Jmax

fitarhJ<-function(dat){
  
  try(Vr<-nls(fVJ.a, start = list(k25=100, Ea=40), data = dat))
  Vr2<-summary(Vr)
  res<-Vr2$coefficients[1:4]
  names(res)[1:4]<-c("Jmax25","EaJ","Jmax25.se","EaJ.se")
  
  #tt<-seq(min(dat$TsK),max(dat$TsK),length=51)
  #predicts<-predictNLS(Vr,newdata=data.frame(TsK=tt),interval="confidence",level=0.95)
  #predicts.df<-data.frame(predicts$summary)
  #predicts.df$TsK<-tt
  
  r<-cor(fitted(Vr),dat$Jmax)
  R2<-r*r
  
  #test for normality of residuals
  rest<-residuals(Vr)
  norm<-shapiro.test(rest)
  S<-norm$statistic
  pvalue<-norm$p.value
  
  param<-c(res,R2)
  names(param)[1:5]<-c("Jmax25","EaJ","Jmax25.se","EaJ.se","R2_Arr")
  return(param)
  
  
  #return(res)
  #return(list(res,predicts.df))
}



#----------------------------------------------------------------------------------------------------------------------
# from Jhon's codes
#- function to fit the June et al. (2004) FPB model for the temperature response of photosynthesis.
#- accepts a dataframe, returns a list with [1] named vector of parameter estiamtes and their se's,
#-   and [2] a dataframe with the predictions and 95% confidence intervals.
fitAvT <- function(dat){
  try(A_Topt <- nls(Photo~ Jref*exp(-1*((Tleaf-Topt)/theta)^2),data=dat,start=list(Jref=20,Topt=25,theta=20)))
  A_Topt2 <- summary(A_Topt)
  results <- A_Topt2$coefficients[1:6]
  names(results)[1:6] <- c("Aopt","Topt","theta","Aopt.se","Topt.se","theta.se")
  
  #TT <- seq(min(dat$Tleaf),max(dat$Tleaf),length=50)
 # predicts <- predictNLS(A_Topt, newdata=data.frame(Tleaf = TT),interval="confidence",level=0.95)
  #predicts.df <- data.frame(predicts$summary)
  #predicts.df$Tleaf <- TT
  
  return(list(results))
  #return(list(results,predicts.df))
}



#-----------------------------------------------------------------------------------------
#function to fit the Quadratic model for the temperature response of photosynthesis.
#model, Anet=Aopt-b(Tleaf-Topt)^2
#- accepts a dataframe, returns a list with [1] named vector of parameter estiamtes and their se's,
#-   and [2] a dataframe with the predictions and 95% confidence intervals.

fitquad<-function(dat){
  try(An<-nls(Photo~Aopt-(b*(Tleaf-Topt)^2),data=dat,start=list(Aopt=max(dat$Photo),Topt=25,b=0.05)))
  A.1<-summary(An)
  results<-(A.1$coefficients[1:6])
  names(results)[1:6] <- c("aopt","topt","b","aopt.se","topt.se","b.se")
  
  #TT.i <- seq(min(dat$Tleaf),max(dat$Tleaf),length=51)
  #predicts <- predictNLS(An, newdata=data.frame(Tleaf = TT.i),interval="confidence",level=0.95)
  #predicts.df <- data.frame(predicts$summary)
  #predicts.df$Tleaf <- TT.i
  
  r<-cor(fitted(An),dat$Photo)
  r2<-r*r
  
  #test for normality of residuals
  rest<-residuals(An)
  norm<-shapiro.test(rest)
  s<-norm$statistic
  pvalue<-norm$p.value
  
  return(c(results,r2,s,pvalue))
  
}

#------------------------------------------------------------------------------------------
#function to fit the Quadratic model for the temperature response of photosynthesis at CI=300.
#model, Anet=Aopt-b(Tleaf-Topt)^2
#- accepts a dataframe, returns a list with [1] named vector of parameter estiamtes and their se's,
#-   and [2] a dataframe with the predictions and 95% confidence intervals.

fitquad300<-function(dat){
  
  try(An<-nls(lowAs~Aopt-(b*(Ts-Topt)^2),data=dat,start=list(Aopt=max(dat$lowAs),Topt=25,b=0.05)))
  
  
  A.1<-summary(An)
  results<-(A.1$coefficients[1:6])
  names(results)[1:6] <- c("aopt","topt","b","aopt.se","topt.se","b.se")
  
  #TT.i <- seq(min(dat$Ts),max(dat$Ts),length=51)
  #predicts <- predictNLS(An, newdata=data.frame(Ts = TT.i),interval="confidence",level=0.95)
  #predicts.df <- data.frame(predicts$summary)
  #predicts.df$Ts <- TT.i
  
  #return(list(results,predicts.df))
  r<-cor(fitted(An),dat$Photo)
  r2<-r*r
  
  #test for normality of residuals
  rest<-residuals(An)
  norm<-shapiro.test(rest)
  s<-norm$statistic
  pvalue<-norm$p.value
  

  return(c(results,r2,s,pvalue))
  
  
  
  
  #return(results)
}

#------------------------------------------------------------------------------------------------------------------------
#function to fit the Quadratic model for the temperature response of photosynthesis at CI=800.
#model, Anet=Aopt-b(Tleaf-Topt)^2
#- accepts a dataframe, returns a list with [1] named vector of parameter estiamtes and their se's,
#-   and [2] a dataframe with the predictions and 95% confidence intervals.


fitquad800<-function(dat){
  
  try(An<-nls(highAs~Aopt-(b*(Ts-Topt)^2),data=dat,start=list(Aopt=max(dat$highAs),Topt=25,b=0.05)))
  
  
  A.1<-summary(An)
  results<-(A.1$coefficients[1:6])
  names(results)[1:6] <- c("aopt","topt","b","aopt.se","topt.se","b.se")
  
  #TT.i <- seq(min(dat$Ts),max(dat$Ts),length=51)
  #predicts <- predictNLS(An, newdata=data.frame(Ts = TT.i),interval="confidence",level=0.95)
  #predicts.df <- data.frame(predicts$summary)
  #predicts.df$Ts <- TT.i
  
  #return(list(results,predicts.df))
  
  #return(list(results,predicts.df))
  r<-cor(fitted(An),dat$highAs)
  r2<-r*r
  
  #test for normality of residuals
  rest<-residuals(An)
  norm<-shapiro.test(rest)
  s<-norm$statistic
  pvalue<-norm$p.value
  
  return(c(results,r2,s,pvalue))
  
  
  
  
  
  #return(results)
}


#---------------------------------------------------------------------------------------------------------------------

#- function to get the confidence intervals of an NLS fit
#https://quantitativeconservationbiology.wordpress.com/2013/07/02/confidence-interval-for-a-model-fitted-with-nls-in-r/

as.lm.nls <- function(object, ...) {
  if (!inherits(object, "nls")) {
    w <- paste("expected object of class nls but got object of class:", 
               paste(class(object), collapse = " "))
    warning(w)
  }
  
  gradient <- object$m$gradient()
  if (is.null(colnames(gradient))) {
    colnames(gradient) <- names(object$m$getPars())
  }
  
  response.name <- if (length(formula(object)) == 2) "0" else 
    as.character(formula(object)[[2]])
  
  lhs <- object$m$lhs()
  L <- data.frame(lhs, gradient)
  names(L)[1] <- response.name
  
  fo <- sprintf("%s ~ %s - 1", response.name, 
                paste(colnames(gradient), collapse = "+"))
  fo <- as.formula(fo, env = as.proto.list(L))
  
  do.call("lm", list(fo, offset = substitute(fitted(object))))
  
}

#----------------------------------------------------------------------------------------------------------------
# Adds error bars to a plot
# From John's codes

adderrorbars <- function(x,y,SE,direction,barlen=0.04,...){
  
  if(length(direction)>1)stop("direction must be of length one.")
  #if(direction == "updown")
  #  direction <- c("up","down")
  if(direction == "rightleft" | direction == "leftright")direction <- c("left","right")
  
  if("up" %in% direction)
    arrows(x0=x, x1=x, y0=y, y1=y+SE, code=3, angle=90, length=barlen,...)
  if("down" %in% direction) 
    arrows(x0=x, x1=x, y0=y, y1=y-SE, code=3, angle=90, length=barlen,...)
  if("updown" %in% direction) 
    arrows(x0=x, x1=x, y0=y+SE, y1=y-SE, code=3, angle=90, length=barlen,...)
  if("left" %in% direction) 
    arrows(x0=x, x1=x-SE, y0=y, y1=y, code=3, angle=90, length=barlen,...)
  if("right" %in% direction)
    arrows(x0=x, x1=x+SE, y0=y, y1=y, code=3, angle=90, length=barlen,...)  
  
}
#----------------------------------------------------------------------------------------------------------------


#is.even <- function(x){ x %% 2 == 0 } 
#aci$c<-ifelse(is.even(aci$Chamber),"Elevated","Ambient")
#----------------------------------------------------------------------------------------------------------------

#function to fit BBOpti stomatal conductance model and to estimate 95% CI for g1


bbopt<-function(gsdata){
  
  fit.nls <- nls(Cond ~ 1.6 * (1 + g1/sqrt(VpdL)) * (Photo/CO2S), 
                 start = list(g1 = 4),data=gsdata)
  
  #g.boot<-nlsBoot(fit.nls,niter=50)
  g <- summary(fit.nls)
  results <- g$coefficients[1:2]
  names(results)[1:2] <- c("g1","Std.err")
  
  
  
  return(results)
}

#----------------------------------------------------------------------------------------------------------------
#function to fit BBOptiFull stomatal conductance model and to estimate 95% CI for g1

bboptful<-function(dat){
  
  fit <- nls(Cond ~ 1.6 * (1 + g1/VpdL^gk) * (Photo/CO2S), 
             start = list(g1 = 4, gk = 0.5),data=dat)
  
  fit2 <- summary(fit)
  ci_g1<-confint2(fit)
  results <- fit2$coefficients[1:4]
  names(results)[1:4] <- c("g1","gk","g1.se","gk.se")
  
  
  #g.boot<-nlsBoot(fit,niter=999)
  
  #vpd <- seq(min(dat$VpdL),max(dat$VpdL),length=50)
  #photo <- seq(min(dat$Photo),max(dat$Photo),length=50)
  #co2s <- seq(min(dat$CO2S),max(dat$CO2S),length=50)
  #predicts <- predictNLS(fit, newdata=data.frame(VpdL = vpd,Photo=photo,CO2S=co2s),interval="confidence",level=0.95)
  #predicts.df <- data.frame(predicts$summary)
  
  
  return(list(results,ci_g1))
}

#------------------------------------------------------------------------------------------------------------------
#function to fit quadratic form to Photosynthetic temperature response
#using dataframe from summaryBy function
#quadratic
fitquadnew<-function(dat){
  try(An<-nls(Photo.mean~Aopt-(b*(Tleaf.mean-Topt)^2),data=dat,start=list(Aopt=max(dat$Photo.mean),Topt=25,b=0.05)))
  A.1<-summary(An)
  results<-(A.1$coefficients[1:6])
  names(results)[1:6] <- c("Aopt","Topt","b","Aopt.se","Topt.se","b.se")
  
  TT.i <- seq(min(dat$Tleaf.mean),max(dat$Tleaf.mean),length=51)
  predicts <- predictNLS(An, newdata=data.frame(Tleaf.mean = TT.i),interval="confidence",level=0.95)
  predicts.df <- data.frame(predicts$summary)
  predicts.df$Tleaf <- TT.i
  
  resids<-residuals(An)
  return(list(results,predicts.df,resids))
  
}
#---------------------------------------------------------------------------------------------------------------------
#June's 2004 model

fitAvTnew <- function(dat){
  try(A_Topt <- nls(Photo.mean~ Jref*exp(-1*((Tleaf.mean-Topt)/theta)^2),data=dat,start=list(Jref=20,Topt=25,theta=20)))
  A_Topt2 <- summary(A_Topt)
  results <- A_Topt2$coefficients[1:6]
  names(results)[1:6] <- c("Aopt","Topt","theta","Aopt.se","Topt.se","theta.se")
  
  TT <- seq(min(dat$Tleaf.mean),max(dat$Tleaf.mean),length=50)
  predicts <- predictNLS(A_Topt, newdata=data.frame(Tleaf.mean = TT),interval="confidence",level=0.95)
  predicts.df <- data.frame(predicts$summary)
  predicts.df$Tleaf <- TT
  
  return(list(results,predicts.df))
}

#-----------------------------------------------------------------------------------------------------------------------
#ITD curve 

fititd<-function(dat){
  try(An<-nls(Photo.mean~(b*(Tleaf.mean-Tmin)*((1-exp(c*(Tleaf.mean-Tmax))))),data=dat,start=list(b=2,Tmin=0,Tmax=50,c=0.03)))
A.1<-summary(An)
results<-(A.1$coefficients[1:6])
  #names(results)[1:6] <- c("Aopt","Topt","b","Aopt.se","Topt.se","b.se")
  
TT.i <- seq(min(dat$Tleaf.mean),max(dat$Tleaf.mean),length=51)
predicts <- predictNLS(An, newdata=data.frame(Tleaf.mean = TT.i),interval="confidence",level=0.95)
  predicts.df <- data.frame(predicts$summary)
  predicts.df$Tleaf <- TT.i
  
  return(list(results,predicts.df))
  
}

#--------------------------------------------------------------------------------------------------------------------------

#functions to fit ACi curves

getr2 <- function(x){
  lmfit <- lm(Ameas ~ Amodel, data=x$df)
  summary(lmfit)$r.squared
}

# function to calculate Topt from fitted quadratic
Topt <- function(pars) {
  -pars[2]/2/pars[3]
}

makecurves <- function(dir, fname) {
  
  data <- read.csv(paste0(dir,"/",fname) )
  
  fits1 <- fitacis(data, "Curve", fitmethod="bilinear", Tcorrect=FALSE,useRd=TRUE)
  fits2 <- fitacis(data, "Curve", fitTPU=T, Tcorrect=FALSE)
  
  # pdf(paste0(dir,"/TPU_ornot.pdf"), width=9, height=5)
  # par(mfrow=c(1,2))
  # for(i in 1:length(fits1)){
  #   plot(fits1[[i]], main=paste(names(fits1)[i], "No TPU", sep=" - "))
  #   plot(fits2[[i]], main=paste(names(fits1)[i], "TPU", sep=" - "))
  # }
  # dev.off()
  
  return(fits2)
}


makedata<- function(path, fname, fit) {
  
  ret <- coef(fit) 
  ret$Jmax[ret$Jmax > 1000] <- NA
  ret$rmse <- sapply(fit, "[[", "RMSE")
  ret$R2 <- sapply(fit, getr2)
  
  # extract fitted values at 275 and 800 ppm Ci
  
  #ret$Ci <- sapply(fit,function(x)min(x$df$Ci))
  ret$resp <- sapply(fit, function(x)x$Photosyn(Ci=50)$ALEAF)
  ret$lowAs <- sapply(fit, function(x)x$Photosyn(Ci=275)$ALEAF)
  ret$highAs <- sapply(fit, function(x)x$Photosyn(Ci=800)$ALEAF)
  ret$elevated <- sapply(fit, function(x)x$Photosyn(Ci=550)$ALEAF)
  
  
  ret$Ci_t <- sapply(fit,function(x)x$Ci_transition)
  
  ret$Ac <- sapply(fit, function(x)x$Photosyn(Ci=275)$Ac)
  ret$Aj <- sapply(fit, function(x)x$Photosyn(Ci=275)$Aj)
  ret$Ap <- sapply(fit, function(x)x$Photosyn(Ci=275)$Ap)
  ret$Ts <- sapply(fit, function(x)mean(x$df$Tleaf))
  ret$JVr<-with(ret,Jmax/Vcmax)
  
  #ret$Condt<-sapply(fit, function(x)mean(x$df$Cond))
  ret$TsK <- ret$Ts+273.15
  ret$Tsq <- ret$Ts * ret$Ts
  ret$maxCi <- sapply(fit,function(x)max(x$df$Ci))
  ret$minCi <- sapply(fit,function(x)min(x$df$Ci))
  
  #to get the limiting process
  ret$lim<-apply(ret[,17:19],1,min)
  
  ret$lim_step<-NA
  ret$lim_step[which(ret$lim==ret$Ac)]<-"C_lim"
  ret$lim_step[which(ret$lim==ret$Aj)]<-"R_lim"
  ret$lim_step[which(ret$lim==ret$Ap)]<-"P_lim"
  ret$lim_step[which(ret$lim==ret$Ac & ret$lim ==ret$Aj)]<-"Co_lim_Ac_Aj"
  
  
  data <- read.csv(paste0(path,"/",fname) )
  ret <- merge(ret,data,by="Curve")
  ret <- ret[!duplicated(ret[,c("Curve")]),]
  ret <- subset(ret,R2 > 0.99)
  
  return(ret)
}


makeplots <- function(dir, data, colby) {
  
  pdf(paste0(dir,"/T_relations.pdf"))
  with(data,plot(Ts,Vcmax,col=colby,main="Vcmax"))
  with(data,plot(Ts,Jmax,col=colby,main="Jmax"))
  with(data,plot(Ts,TPU,col=colby,main="TPU"))
  with(data,plot(Ts,Rd,col=colby,main="Rday"))
  with(data,plot(Ts,resp,col=colby,main="A at Ci = 50 ppm"))
  with(data,plot(Ts,lowAs,col=colby,main="A at Ci = 300 ppm"))
  with(data,plot(Ts,highAs,col=colby,main="A at Ci = 800 ppm"))
  dev.off()
  
}

#in this function, R2 for selecting curves was reduced to 0.9.
#this is only used for SAVANNA species and only 2 curves have R2 <0.99 (both > 0.9)
#these two curves measured at 25C and only these two available for T=25C

makedata.1<- function(path, fname, fit) {
  
  ret <- coef(fit) 
  ret$Jmax[ret$Jmax > 1000] <- NA
  ret$rmse <- sapply(fit, "[[", "RMSE")
  ret$R2 <- sapply(fit, getr2)
  
  # extract fitted values at 275 and 800 ppm Ci
  
  #ret$Ci <- sapply(fit,function(x)min(x$df$Ci))
  ret$resp <- sapply(fit, function(x)x$Photosyn(Ci=50)$ALEAF)
  ret$lowAs <- sapply(fit, function(x)x$Photosyn(Ci=275)$ALEAF)
  ret$highAs <- sapply(fit, function(x)x$Photosyn(Ci=800)$ALEAF)
  ret$elevated <- sapply(fit, function(x)x$Photosyn(Ci=550)$ALEAF)
  
  
  ret$Ci_t <- sapply(fit,function(x)x$Ci_transition)
  
  ret$Ac <- sapply(fit, function(x)x$Photosyn(Ci=275)$Ac)
  ret$Aj <- sapply(fit, function(x)x$Photosyn(Ci=275)$Aj)
  ret$Ap <- sapply(fit, function(x)x$Photosyn(Ci=275)$Ap)
  ret$Ts <- sapply(fit, function(x)mean(x$df$Tleaf))
  ret$JVr<-with(ret,Jmax/Vcmax)
  
  #ret$Condt<-sapply(fit, function(x)mean(x$df$Cond))
  ret$TsK <- ret$Ts+273.15
  ret$Tsq <- ret$Ts * ret$Ts
  ret$maxCi <- sapply(fit,function(x)max(x$df$Ci))
  ret$minCi <- sapply(fit,function(x)min(x$df$Ci))
  
  #to get the limiting process
  ret$lim<-apply(ret[,17:19],1,min)
  
  ret$lim_step<-NA
  ret$lim_step[which(ret$lim==ret$Ac)]<-"C_lim"
  ret$lim_step[which(ret$lim==ret$Aj)]<-"R_lim"
  ret$lim_step[which(ret$lim==ret$Ap)]<-"P_lim"
  ret$lim_step[which(ret$lim==ret$Ac & ret$lim ==ret$Aj)]<-"Co_lim_Ac_Aj"
  
  
  data <- read.csv(paste0(path,"/",fname) )
  ret <- merge(ret,data,by="Curve")
  ret <- ret[!duplicated(ret[,c("Curve")]),]
  ret <- subset(ret,R2 > 0.9)
  
  return(ret)
}



#----------------------------------------------------------------------------------------------------------------------

#function to fit linear mix models to estimate Topt for photosynthesis from quadratic equation

#---------------------------------------------------------------------------------------------------------------------
#dat=dataframe
#As=dependent variable (Photo, lowAs, highAs.....)
#rand=random effect

fit.boot<-function(dat,As,rand){
  
  q<-lmer(As~Ts+Tsq+(1|rand),data=dat)
  
  #function to get Topt for photosynthesis and Rate at Topt
  
  topt<-function(q){
    q.1<-summary(q)
    results<-(q.1$coefficients[1:3])
    Topt<--results[2]/2/results[3]
    Aopt<-results[1]+results[2]*Topt+results[3]*Topt^2
    
    
    return(c(Topt,Aopt))
  }
  
  #to get bootstrap SE
  bootfit2 <- bootMer(q, FUN=topt, re.form=NA,
                      nsim=999)
  
  #to get bootstrap CI
  
  lci <- apply(bootfit2$t, 2, quantile, 0.025)
  uci <- apply(bootfit2$t, 2, quantile, 0.975)
  
  return(list(bootfit2,lci,uci))
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

fititd<-function(dat){
  try(An<-nls(Photo~(b*(Tleaf-Tmin)*((1-exp(c*(Tleaf-Tmax))))),data=dat,start=list(b=2,Tmin=0,Tmax=40,c=0.03)))
  A.1<-summary(An)
  results<-(A.1$coefficients[1:6])
  #names(results)[1:6] <- c("Aopt","Topt","b","Aopt.se","Topt.se","b.se")
  return(results)
  
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

acifit <- function(path, fname) {
  
  data <- read.csv(paste0(path,"/",fname) )
  
  fits1 <- fitacis(data, "Curve", fitmethod="bilinear", Tcorrect=FALSE)
  fits2 <- fitacis(data, "Curve", fitTPU=TRUE, Tcorrect=FALSE)
  
  
  return(fits2)
}

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


#function to capital first letter in a character string

capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}



#get mean Rday and Photosynthesis @ Ci=275 with one factor variable
get_summary_one<-function(data,from,to,by1){
  
  data$by1 <- data[[by1]]
  
  data.1<-data[data$Ts >= from & data$Ts <= to,]
  
  #data$TargetTemp<-cut(data$Ts,breaks=breaks,
  #labels=FALSE)
  Rday<-summaryBy(JVr+lowAs+Rd+Ts~by1,data=data.1,FUN=c(mean,std.error),na.rm=T)
  
  Rday.1<-Rday[-c(5,9)]
  names(Rday.1)[1:7]<-c(by1,"JVr","Photo_275","Rday","JVr_SE","Photo_275_SE","Rday_SE")
  return(Rday.1)
}
#get_summary_one(damz,by1="Season",from=23,to=27)

#get mean Rday and Photosynthesis @ Ci=275 with two factor variable

get_summary_two<-function(data,by1,by2,from,to){
  
  data$by1 <- data[[by1]]
  data$by2 <- data[[by2]]
  
  data.1<-data[data$Ts >= from & data$Ts <= to,]
  Rday<-summaryBy(JVr+lowAs+Rd+Ts~by1+by2,data=data.1,FUN=c(mean,std.error),na.rm=TRUE)
  
  Rday.2<-Rday[-c(6,10)]
  names(Rday.2)[1:8]<-c(by1,by2,"JVr","Photo_275","Rday","JVr_SE","Photo_275_SE","Rday_SE")
  return(Rday.2)
}



get_topts<-function(tfits){
  fc<-list()
  for(i in 1:length(tfits)){
    fc[[i]] <- tfits[[i]]
    topts<-do.call(rbind,fc)
    topts<-data.frame(topts)
    
    
  }
  return(topts)
}

#functions to pullout names of lists

get_names<-function(list.1){
  Species<-list()
  Season<-list()
  
  Treatment<-strsplit(names(list.1)," ")
  
  for(i in 1:length(list.1)){
    Species[i]<-Treatment[[i]][[1]]
    Season[i]<-Treatment[[i]][[2]] 
    
    #Species<-do.call(cbind,(spp))
    #Season<-do.call(cbind,(sea))
    
    topts<-data.frame((cbind(Species,Season)))
    topts.1<-cbind(ldply(topts$Species),ldply(topts$Season))
  }
  return(topts.1)
}

#function to fit temperature response of photosynthesis (mix model)

fit.nlme<-function(dat,random,yvar){
  
  dat$rand <- dat[[random]]
  dat$yvar <- dat[[yvar]]
  
  testfit2<-nlme(yvar~Aopt-(b*(Tleaf-Topt)^2),fixed=list(Aopt + Topt + b ~ 1),random = Aopt+Topt ~ 1 | rand,
                 start=list(fixed=c(Aopt=max(dat$Photo),Topt=25,b=0.05)),data=dat)
  
  test<-testfit2$tTable
  
  aopt<-summary(testfit2)$tTable[[1]]
  topt<-summary(testfit2)$tTable[[2]]
  b<-summary(testfit2)$tTable[[3]]
  aopt.se<-summary(testfit2)$tTable[[4]]
  topt.se<-summary(testfit2)$tTable[[5]]
  b.se<-summary(testfit2)$tTable[[6]]
  
  #to get R2 between fitted and observed photosynthesis
  r<-cor(fitted(testfit2),dat$Photo)
  r2<-r*r
  
  #test for normality of residuals
  rest<-residuals(testfit2)
  norm<-shapiro.test(rest)
  s<-norm$statistic
  pvalue<-norm$p.value
  
  param<-cbind(aopt,topt,b,aopt.se,topt.se,b.se)
  
  names(param)[1:6] <- c("aopt","topt","b","aopt.se","topt.se","b.se")
  
  #return(param)
  
  
  return(c(param,r2,s,pvalue))
  
}

#read and process Mike Aspinwall's data

#corm<-read.csv(paste0(path,"/Data/mike_aspinwall_corymbia_calophylla.csv"))

#corm$Curve <- c(1)
#Curve <- c()
#count <- 1  
#for (i in 2:length(corm$CO2R)){
  #ifelse(abs(corm$CO2S[i-1] - corm$CO2S[i])<500,count <- count,count <- count + 1)
 # Curve[i] <- count 
#}
#corm$Curve[2:length(corm$CO2S)] <- na.omit(Curve)

#corm<-subset(corm,!is.na(Photo))
#write.csv(corm,paste0(path,"/Data/mike_aspinwall_corymbia_calophylla.V1.csv"),row.names=FALSE,sep=",")


#functions to fit temperature response of Vcmax and Jmax in mix model framework
#These functions use extimated Vcmax and Jmax returned from the plantecophys package.
#fits peaked model and returen temperature response parameters, R2 (squared correlation coefficient between fitted and 
#observed Vcmax or Jmax and test for normality of residuals (Sharpiro-Whilks test))

fitvcmax_mm<-function(dat,random,return=c("Arr","Peak"),start=list(fixed=c(k25=100, Ea=60, delS = 0.64))){
  
  dat$rand <- dat[[random]]
  return <- match.arg(return)
  
  try(Vc<-nlme(fVc, fixed=list(k25 + Ea + delS ~ 1),random = k25 + Ea ~ 1 | rand,
               start = start, data = dat))
  
  
  try(Vr<-nlme(fVc.a, fixed=list(k25 + Ea ~ 1),random = k25 + Ea ~ 1 | rand,
               start = list(fixed=c(k25=100, Ea=60)), data = dat))
  
  #test for best fitted function
  
  an<-anova(Vr,Vc)
  AIC_Arr<-an$AIC[1]
  AIC_Peak<-an$AIC[2]
  prob_V<-an[[9]][[2]]
  
  Vcmax25_a<-summary(Vc)$tTable[[1]]
  Ea_a<-summary(Vc)$tTable[[2]]
  delS<-summary(Vc)$tTable[[3]]
  Vcmax25.se_a<-summary(Vc)$tTable[[4]]
  Ea.se_a<-summary(Vc)$tTable[[5]]
  delS.se<-summary(Vc)$tTable[[6]]
  
  para_peak<-data.frame(cbind(Vcmax25_a,Ea_a,delS,Vcmax25.se_a,Ea.se_a,delS.se))
  names(para_peak)[1:6]<-c("Vcmax25","EaV","delsV","Vcmax25.se","EaV.se","delsV.se")
  
  #pull out fitted parameters of Arrh model
  Vcmax25_b<-summary(Vr)$tTable[[1]]
  Ea_b<-summary(Vr)$tTable[[2]]
  Vcmax25.se_b<-summary(Vr)$tTable[[3]]
  Ea.se_b<-summary(Vr)$tTable[[4]]
  para_arr<-data.frame(cbind(Vcmax25_b,Ea_b,Vcmax25.se_b,Ea.se_b))
  names(para_arr)[1:4]<-c("Vcmax25","EaV","Vcmax25.se","EaV.se")
  
  
  #R2 
  r1<-cor(fitted(Vc),dat$Vcmax)
  R2_Peak<-r1*r1
  
  r2<-cor(fitted(Vr),dat$Vcmax)
  R2_Arr<-r2*r2
  
  #test for normality of residuals
  rest<-residuals(Vc)
  norm<-shapiro.test(rest)
  S<-norm$statistic
  pvalue<-norm$p.value
  
  #to calculate Topt for Vcmax
  topt<-200/(delS-0.008314*log(Ea_a/(200-Ea_a)))
  ToptV<-topt-273.15
  
  param_peaked<-cbind(para_peak,ToptV,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_V)
  param_arr<-cbind(para_arr,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_V)
  
  if(return == "Peak")return(param_peaked)
  if(return == "Arr")return(param_arr)
  
}

# function to fit temperature response of Vcmax (without random effects)

fitpeaked<-function(dat,return=c("Arr","Peak"),start=list(k25=100, Ea=60, delS = 0.64)){
  
  return <- match.arg(return)
  try(Vc<-nls(fVc, start=start, data=dat))
  Vc2<-summary(Vc)
  res1<-Vc2$coefficients[1:6]
  names(res1)[1:6]<-c("Vcmax25","EaV","delsV","Vcmax25.se","EaV.se","delsV.se")
  
  try(Vr<-nls(fVc.a, start = list(k25=100, Ea=40), data = dat))
  Vr2<-summary(Vr)
  res<-Vr2$coefficients[1:4]
  names(res)[1:4]<-c("Vcmax25","EaV","Vcmax25.se","EaV.se")
  
  
  an<-anova(Vr,Vc)
  AIC_Arr<-AIC(Vr)
  AIC_Peak<-AIC(Vc)
  prob_V<-an[[6]][[2]]
  
  
  r1<-cor(fitted(Vc),dat$Vcmax)
  R2_Peak<-r1*r1
  
  r2<-cor(fitted(Vr),dat$Vcmax)
  R2_Arr<-r2*r2
  
  #test for normality of residuals
  rest<-residuals(Vc)
  norm<-shapiro.test(rest)
  s<-norm$statistic
  pvalue<-norm$p.value
  
  topt<-200/(res1[[3]]-0.008314*log(res1[[2]]/(200-res1[[2]])))
  Topt<-topt-273.15
  
  param_peak<-c(res1,Topt,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_V)
  names(param_peak)[1:12]<-c("Vcmax25","EaV","delsV","Vcmax25.se","EaV.se","delsV.se","ToptV","R2_Arr",
                             "R2_Peak","AIC_Arr","AIC_Peak","prob_V")
  
  param_arr<-c(res,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_V)
  names(param_arr)[1:9]<-c("Vcmax25","EaV","Vcmax25.se","EaV.se","R2_Arr",
                           "R2_Peak","AIC_Arr","AIC_Peak","prob_V")
  
  if(return == "Peak")return(param_peak)
  if(return == "Arr")return(param_arr)
  
  
}

# function to fit temperature response of jmax (mix model)

fitjmax_mm<-function(dat,random,return=c("Arr","Peak"),start = list(fixed=c(k25=100, Ea=60, delS = 0.64))){
  return <- match.arg(return)
  dat$rand <- dat[[random]]
  
  try(Vc<-nlme(fVJ, fixed=list(k25 + Ea + delS ~ 1),random = k25 + Ea ~ 1 | rand,
               start = start, data = dat))
  
  try(Vr<-nlme(fVJ.a, fixed=list(k25 + Ea ~ 1),random = k25 + Ea ~ 1 | rand,
               start = list(fixed=c(k25=100, Ea=60)), data = dat))
  
  
  an<-anova(Vr,Vc)
  AIC_Arr<-an$AIC[1]
  AIC_Peak<-an$AIC[2]
  prob_J<-an[[9]][[2]]
  
  
  Jmax25<-summary(Vc)$tTable[[1]]
  EaJ<-summary(Vc)$tTable[[2]]
  delsJ<-summary(Vc)$tTable[[3]]
  Jmax25.se<-summary(Vc)$tTable[[4]]
  EaJ.se<-summary(Vc)$tTable[[5]]
  delsJ.se<-summary(Vc)$tTable[[6]]
  
  
  Jmax25_b<-summary(Vr)$tTable[[1]]
  Ea_b<-summary(Vr)$tTable[[2]]
  Jmax25.se_b<-summary(Vr)$tTable[[4]]
  Ea.se_b<-summary(Vr)$tTable[[5]]
  para_arr<-data.frame(cbind(Jmax25_b,Ea_b,Jmax25.se_b,Ea.se_b))
  names(para_arr)[1:4]<-c("Jmax25","EaJ","Jmax25.se","EaJ.se")
  
  
  #R2 
  r1<-cor(fitted(Vc),dat$Vcmax)
  R2_Peak<-r1*r1
  
  r2<-cor(fitted(Vr),dat$Vcmax)
  R2_Arr<-r2*r2
  
  #test for normality of residuals
  rest<-residuals(Vc)
  norm<-shapiro.test(rest)
  S<-norm$statistic
  pvalue<-norm$p.value
  
  topt<-200/(delsJ-0.008314*log(EaJ/(200-EaJ)))
  ToptJ<-topt-273.15
  
  
  param_peak<-cbind(Jmax25,EaJ,delsJ,Jmax25.se,EaJ.se,delsJ.se,ToptJ,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_J)
  param_arr<-cbind(para_arr,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_J)
  
  if(return == "Peak")return(param_peak)
  if(return == "Arr")return(param_arr)
  
}


# function to fit temperature response of jmax (without random effects)

fitpeakedJ<-function(dat,return=c("Arr","Peak"),start=list(k25=100, Ea=60, delS = 0.64)){
  return <- match.arg(return)
  
  try(Vj<-nls(fVJ, start = start, data = dat))
  Vj2<-summary(Vj)
  res1<-Vj2$coefficients[1:6]
  names(res1)[1:6]<-c("Jmax25","Ea","delsJ","Jmax25.se","EaJ.se","delsJ.se")
  
  try(Vr<-nls(fVJ.a, start = list(k25=100, Ea=40), data = dat))
  Vr2<-summary(Vr)
  res<-Vr2$coefficients[1:4]
  names(res)[1:4]<-c("Vcmax25","EaV","Vcmax25.se","EaV.se")
  
  
  an<-anova(Vr,Vj)
  AIC_Arr<-AIC(Vr)
  AIC_Peak<-AIC(Vj)
  prob_J<-an[[6]][[2]]
  
  
  r1<-cor(fitted(Vj),dat$Jmax)
  R2_Peak<-r1*r1
  
  r2<-cor(fitted(Vr),dat$Jmax)
  R2_Arr<-r2*r2
  
  
  
  #test for normality of residuals
  rest<-residuals(Vj)
  norm<-shapiro.test(rest)
  s<-norm$statistic
  pvalue<-norm$p.value
  
  topt<-200/(res1[[3]]-0.008314*log(res1[[2]]/(200-res1[[2]])))
  Topt<-topt-273.15
  
  
  param_peak<-c(res1,Topt,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_J)
  names(param_peak)[1:12]<-c("Jmax25","EaJ","delsJ","Jmax25.se","EaJ.se","delsJ.se","ToptJ","R2_Arr","R2_Peak",
                             "AIC_Arr","AIC_Peak","prob_J")
  
  param_arr<-c(res,R2_Arr,R2_Peak,AIC_Arr,AIC_Peak,prob_J)
  names(param_arr)[1:9]<-c("Jmax25","EaJ","Jmax25.se","EaJ.se","R2_Arr",
                           "R2_Peak","AIC_Arr","AIC_Peak","prob_J")
  
  if(return == "Peak")return(param_peak)
  if(return == "Arr")return(param_arr)
  
  #return(c(res,r2,s,pvalue))
  #return(list(res))
}

#read and process Amazon gs data

read_and_convert_licor <- function(path){
  
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  files.summer <- list.files(pattern="\\.csv")
  allxlfiles.s<-list()
  
  for (i in seq_along(files.summer)){
    
    
    dat <- suppressWarnings(read.csv(files.summer[i]))
    
    dat <- as.data.frame(dat)
    
    dat <- subset(dat,!is.na(as.numeric(dat$Obs)))
    
    allxlfiles.s[[i]] <- dat
    
  }
  
  summer<-do.call(rbind.fill,allxlfiles.s)
  
  summer<-summer[c("Photo","Ci","Tleaf","Cond","CO2S","CO2R","VpdL","Data","PARi")]
  summer<-subset(summer,summer$PARi>990)
  summer<-subset(summer,summer$CO2R>390 & summer$CO2R<425)
  
  #to remove one outliers (negative gs)
  summer<-subset(summer,summer$Cond>0)
  
  
  
  return(summer)
  
}


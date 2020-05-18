
# library("plantecophys")
library("nlstools")
# library(ggplot2)
# library(plover)
# library(RColorBrewer)
# library(plotBy)

# functions

#-----------------------------------------------------------------------------------------
#- function to read and process the Asat and AQ datasets
getAQ <- function(path="Data/Glasshouse_DRAKE_EUTE_THERMAL-NICHE/data"){
  
  #- read in the first set of measurements
  aq1 <-read.csv(paste(path,"GHS39_GREAT_MAIN_GX-AQ_20160202-20160203_L3.csv",sep="/"))
  aq1$campaign = 1
  
  #- read in the second set (prov B only!)
  aq2 <-read.csv(paste(path,"GHS39_GREAT_MAIN_GX-AQ_20160225-20160226_L2.csv",sep="/"))
  aq2$campaign = 2
  
  aq <- rbind(aq1,aq2)
  #names(aq)[1:2] <- tolower(names(aq)[1:2])
  #aq$prov <- as.factor(substr(aq$pot,start=1,stop=1))
  #aq$room <- as.factor(aq$room)
  aq$prov_trt <- as.factor(paste(aq$Prov,aq$Room,sep="-"))
  
  #- assign drought treatments
  #aq$Water_trt <- "wet"
  #aq$Water_trt[grep("Bd",aq$pot)] <- "dry"
  aq$W_treatment <- factor(aq$W_treatment,levels=c("w","d"))
  
  #- assign light levels to a factor variable
  aq$LightFac <- NA
  aq$LightFac[which(aq$PARi<150)] <- 1
  aq$LightFac[which(aq$PARi>400 & aq$PARi<600)] <- 2
  aq$LightFac[which(aq$PARi>800 & aq$PARi<1100)] <- 3
  aq$LightFac[which(aq$PARi>1200)] <- 4
  aq$LightFac <- as.factor(aq$LightFac)
  
  #- a few pots have strange data. UWS3 didn't seem to actually drop to low PARi on three cases.
  #-   so remove those bad data. This should remove 3 datapoints!
  aq <- aq[complete.cases(aq),]
  
  #- assign the temperature levels
  aq$TleafFac <- cut(aq$Tleaf,breaks=c(15,22,26,29,34,37,45),labels=1:6)
  
  #- merge in the location factor
  key2 <- data.frame(Prov=as.factor(LETTERS[1:3]),location= factor(c("Cold-edge","Warm-edge","Central"),
                                                                   levels=c("Cold-edge","Central","Warm-edge")))
  aq2 <- merge(aq,key2,by="Prov")
  
  #- average across replicate logs
  aq.means <- summaryBy(.~Room+Pot+Unit+Prov+prov_trt+location+W_treatment+LightFac+TleafFac+campaign,data=aq2,FUN=mean,keep.names=T)
  
  #- remove a crazy datapoint
  torm <- which(aq.means$prov_trt=="A-3" & round(aq.means$CTleaf,2)==28.75)
  aq.means[torm,"Photo"] <- NA 
  
  return(aq.means)
}
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#- function to read and process the temperature response curves of photosynthesis
getAvT <- function(path="Data/Glasshouse_DRAKE_EUTE_THERMAL-NICHE/data"){
  
  avt <-read.csv(paste(path,"GHS39_GREAT_MAIN_GX-AVT_20160205_L3.csv",sep="/"))
  avt$Room <- as.factor(avt$Room)
  avt$prov_trt <- as.factor(paste(avt$Prov,avt$Room,sep="-"))
  
  #- assign drought treatments
  avt$W_treatment <- factor(avt$W_treatment,levels=c("w","d"))
  
  #- assign light levels to a factor variable
  avt$LightFac <- NA
  avt$LightFac[which(avt$PARi<150)] <- 1
  avt$LightFac[which(avt$PARi>1200)] <- 4
  avt$LightFac <- as.factor(avt$LightFac)
  
  #- assign the temperature levels
  avt$TleafFac <- cut(avt$Tleaf,breaks=c(15,22,26,29,34,37,45),labels=1:6)
  
  key2 <- data.frame(Prov=as.factor(LETTERS[1:3]),location= factor(c("Cold-edge","Warm-edge","Central"),
                                                                   levels=c("Cold-edge","Central","Warm-edge")))
  avt2 <- merge(avt,key2,by="Prov")
  
  #avt3 <- summaryBy(.~Room+Code+Unit+Prov+prov_trt+W_treatment+LightFac+TleafFac,data=avt2,FUN=mean,keep.names=T)
  return(avt2)
}
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


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


makedata.1<- function(path, fname, fit) {
  
  ret <- coef(fit) 
  ret$Jmax[ret$Jmax > 1000] <- NA
  ret$rmse <- sapply(fit, "[[", "RMSE")
  ret$R2 <- sapply(fit, getr2)

  #ret$Condt<-sapply(fit, function(x)mean(x$df$Cond))
  ret$Ts <- sapply(fit, function(x)mean(x$df$Tleaf))
  ret$TsK <- ret$Ts+273.15
  
  data <- read.csv(paste0(path,"/",fname) )
  ret <- merge(ret,data,by="Curve")
  ret <- ret[!duplicated(ret[,c("Curve")]),]
  ret <- subset(ret,R2 > 0.90)
  
  return(ret)
}


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


# fit temperature response of Vcmax and get CI of parameters


fVc <- as.formula(Vcmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)) * 
                    (1+exp((298.15*delS - 200)/(298.15*0.008314))) / 
                    (1+exp((TsK*delS-200)/(TsK*0.008314))))

fVJ <- as.formula(Jmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)) * 
                    (1+exp((298.15*delS - 200)/(298.15*0.008314))) / 
                    (1+exp((TsK*delS-200)/(TsK*0.008314))))

fVc.a<-as.formula(Vcmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)))
fVJ.a<-as.formula(Jmax ~ k25 * exp((Ea*(TsK - 298.15))/(298.15*0.008314*TsK)))


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

#------------------------------------------------------------------------------------------------------------------

#function to add a fitted line (linear regression model) to a plot 

fit_line<-function(data,yvar,xvar,byvar=NULL,fitoneline=FALSE,linecol){
  
  data$yvar<-data[[yvar]]
  data$xvar<-data[[xvar]]
  
  
  
  if(!fitoneline){
    data$byvar<-data[[byvar]]
    
    lmlis1 <- lmList(yvar ~ xvar|byvar, data=data)
    liscoef <- coef(lmlis1)
    dfr<-split(data,factor(data[[byvar]]))
    
    for (i in 1:length(dfr)) {
      
      xmin <- min(dfr[[i]][[xvar]])
      xmax <- max(dfr[[i]][[xvar]])
      
      ablineclip(liscoef[i, 1], liscoef[i, 2], x1 = xmin, x2 = xmax,col=linecol,lwd=3)
    }
  }
  
  
  else{
    
    
    lmfit<-lm(yvar~xvar,data=data)
    
    xmin <- min(data[[xvar]])
    xmax <- max(data[[xvar]])
    
    ablineclip(a=coef(lmfit)[[1]],b=coef(lmfit)[[2]],x1=xmin,x2=xmax,col=linecol,lwd=3)
    
  }
  
  
}

#---------------------------------------------------------------------------------------------------------------

# function to fit Medlyn et al stomatal conductance model and get SEs

fitStom<-function(data){
  
  myfit<-fitBB(data,gsmodel="BBOptiFull",fitg0=F)
  coefs<-summary(myfit$fit)$parameters
  g1<-coefs[1]
  g1.se<-coefs[3]
  param<-data.frame(cbind(g1,g1.se))
  return(param)
}

#---------------------------------------------------------------------------------------------------------------

# some other routinely needed functions

#some functions to get  temperature response of gammastar

.Rgas <- function()8.314
Tk <- function(x)x+273.15

# Arrhenius model
arrh <- function(Tleaf,Ea,Tb){
  exp((Ea * (Tk(Tleaf) - Tk(Tb))) / (Tk(Tb) * .Rgas() * Tk(Tleaf))) 
}

# function to get temperature response of GammaStar

TGammaStar <- function(Tleaf, Patm,
                       Egamma=37830.0, 
                       value25=42.75,Tb=25){  
  
  value25*arrh(Tleaf,Egamma,Tb=Tb)*Patm/100
}

# function to get Rday or Rdark from temperature response parameters (Arrhenius model)
TRday <- function(Tleaf, Patm=100,
                  Egamma=46390, 
                  value25=.68,Tb=22){  
  
  value25*arrh(Tleaf,Egamma,Tb)*Patm/100
}


# function to get Rday or Rdark from temperature response parameters (Q10 model)
TRdayQ10<-function(Rref=0.68,Q10=2.1,Tref=22,Tleaf){
  Rday<-Rref*Q10^((Tleaf-Tref)/10)
  
}

#-----------------------------------------------------------------------------------------------------------------

#function to estimate alpha from linear regression coefficients 

get_alpha<-function(data){
  
  lm_a<-lm(Photo~xt,data=data)
  
  alpha<-summary(lm_a)$coefficients[2]
  alpha.se<-summary(lm_a)$coefficients[4]
  param<-data.frame(cbind(alpha,alpha.se))
  return(param)
  
}

#-----------------------------------------------------------------------------------------------------------------

#function to fit temperature response of photosynthesis

fit.nlme<-function(dat,random,yvar,xvar){
  
  dat$rand <- dat[[random]]
  dat$yvar <- dat[[yvar]]
  dat$xvar <- dat[[xvar]]
  
  testfit2<-nlme(yvar~Aopt-(b*(xvar-Topt)^2),fixed=list(Aopt + Topt + b ~ 1),random = Aopt+Topt ~ 1 | rand,
                 start=list(fixed=c(Aopt=max(dat[yvar]),Topt=25,b=0.05)),data=dat)
  
  test<-testfit2$tTable
  
  aopt<-summary(testfit2)$tTable[[1]]
  topt<-summary(testfit2)$tTable[[2]]
  b<-summary(testfit2)$tTable[[3]]
  aopt.se<-summary(testfit2)$tTable[[4]]
  topt.se<-summary(testfit2)$tTable[[5]]
  b.se<-summary(testfit2)$tTable[[6]]
  
  #to get R2 between fitted and observed photosynthesis
  r<-cor(fitted(testfit2),dat[yvar])
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


#-------------------------------------------------------------------------------------------------------------------
#--- estimate root and stem respiration using Q10 model
#--- use measured basal rates at 25 C and Q10 same as leaf dark respiration (Q10=2.1 averaged across provenance)


#--- function to calculate stem and root respiration rates using 15 min temperature measurements
# met<-getRes(met,Rref_S=0.2,Rref_R=0.4,Q10=2.1)
getRes<-function(met,Rref_S,Rref_R,Q10){
  
  met$Rstem<- with(met,Rref_S*Q10^((Tair-25)/10))
  met$Rroot<- with(met,Rref_R*Q10^((Tair-25)/10))
  
  #R_component<-met[c("DateTime_hr", "Bay", "Room","Date","Tair","Rstem","Rroot")]
  return(met)
}


#_-------------------------------------------------------------------------------------------------------------

#-- Function to fit smooth splin to data
interp_var <- function(d, datesout,var){
  
  sp <- spline(x=d$Date, y=d[[var]], xout=datesout)
  out<-data.frame(Date=datesout, Room=unique(d$Room), var=sp$y)
  return(out)
}

#-------------------------------------------------------------------------------------------------------------


#-- function to fit linear model and get parameters
fit_lm<-function(data,yvar,xvar){
  
  
  data$yvar <- data[[yvar]]
  data$xvar <- data[[xvar]]
  
  #mod<-lmer(yvar~xvar+(1|Growth_Condition)+(1|Age_of_Plants),data=data)
  mod<-lm(yvar~xvar,data=data)
  summary(mod)$coefficients
  
  Intercept<-summary(mod)$coefficients[[1]]
  Intercept.se<-summary(mod)$coefficients[[3]]
  Slope<-summary(mod)$coefficients[[2]]
  Slope.se<-summary(mod)$coefficients[[4]]
  Pr_Int<-summary(mod)$coefficients[[7]]
  Pr_Sl<-summary(mod)$coefficients[[8]]
  
  
  Rsq<-summary(mod)$adj.r.squared
  AIC<-AIC(mod)
  Pvalue<-glance(mod)$p.value 
  
  parm<-data.frame(cbind(Intercept,Intercept.se,Slope,Slope.se,Pr_Int, Pr_Sl,Rsq,AIC,Pvalue))
  #names(parm)[1:8]<-c("Intercept","Intercept.se","Slope","Slope.se","Chisq","Pr","Rsq","AIC")
  #parm<-data.frame(parm)
  
  return(parm)
  
}

#-------------------------------------------------------------------------------------------------

# function to copy and paste files


copy_and_paste<-function(from="C:/Repos/DA_GREAT_watered/output",
                         to="C:/Repos/GREATDA/output",lastn="Figure.png"){
  
  list.of.files<-list.files(from, pattern=lastn, all.files=T, full.names=T)
  file.copy(list.of.files, to,overwrite=T)
  
  list.of.files.n<-list.files(to, pattern=lastn, all.files=T, full.names=T)
  
}
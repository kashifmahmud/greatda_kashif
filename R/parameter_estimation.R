#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------

#--- this script estimate the required physiological parameters for model daily carbon gain of Eucalyptus tereticornis prevenances
#--- grown in six growth temperatures. 
#--- parameters : 1. g1 (Medlyn et al stomatal conductance model; both short-term and long-term temperature responces)
#               : 2. alpha (quantum yield of electon transport; by two methods)
#               : 3. Vcmax, Jmax, EaV, EaJ, deslSv, delSj

#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------


# library(devtools)
# install_bitbucket("remkoduursma/plantecophys")
# library(plantecophys)
#-------------------------------------------------------
#- get the photosynthesis vs. temperature data
#- this dataset contains measurements from two campaigns on 2016/02/2 (campaign 1) and 2016/02/25 (campaign 2)
#- in campaign 2, provenance B only measured
#- measurements at different PAR levels 100,500,100,1800 mumols
#avt <- getAvT()

#--- read in-situ photosynthesis measurements (long-term photosynthesis data)

avt.l<-getAQ()

#-------------------------------------------------------

#- get data for campaign 1 

#avt.1<-subset(avt.l,avt.l$campaign==1)


avt.l.list <- subset(avt.l, avt.l$campaign==1 & avt.l$LightFac %in% c(4) & avt.l$W_treatment=="w") 
# avt.l.list <- subset(avt.l, avt.l$LightFac %in% c(4)) 

#-------------------------------------------------------

#fit Medlyn et al 2011 stomatal conductance model
#g0 assumed to be = 0

tofit.l <- split(avt.l.list,avt.l.list$Room)

fitStom<-function(data){
  
  myfit<-fitBB(data,gsmodel="BBOptiFull",fitg0=F)
  coefs<-summary(myfit$fit)$parameters
  g1<-coefs[1]
  g1.se<-coefs[3]
  param<-data.frame(cbind(g1,g1.se))
  return(param)
  }

g1fits<- get_topts(lapply(tofit.l,FUN=fitStom))
g1fits$Room<-c(seq(1,6,1))
g1fits$Tgrowth<-c(18,21.5,25,28.5,32,35.5)

write.csv(g1fits,file = "Parameters/great_g1_long_term.csv",row.names=FALSE,sep=",")
#----------------------------------------------------------------------------------------

#- get data for campaign 2 

# avt.c2<-subset(avt,avt$campaign==2)
# 
# avt.3 <- subset(avt.c2, LightFac %in% c(4) & avt.c2$W_treatment=="w") 
# 
# 
# g1fits<- get_topts(lapply(tofit.2,FUN=fitStom))
# g1fits$Room<-c(seq(1,6,1))
# g1fits$Tgrowth<-c(18,21.5,25,28.5,32,35.5)

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

#- estimate g1 from ACi data

#- read ACi data
path<-getwd()
gr_aci<-read.csv(paste0(path,"/Data/Glasshouse_DRAKE_EUTE_THERMAL-NICHE/data/GHS39_GREAT_MAIN_ACiT_20160216-20160227_L0.csv"))

#- get ambient CO2 levels (first observation)

gr_aci.a<-gr_aci[firstobs(~Curve,data=gr_aci),]

gr_aci.a$gs_pred<-with(gr_aci.a,Photo/(CO2S*sqrt(VpdL)))
gr_aci.a<-subset(gr_aci.a,!is.na(gs_pred) & gr_aci.a$Cond<2) #remove three datapoints with cond>2
gr_aci.a$TleafFac <-cut(gr_aci.a$Tleaf,breaks=c(15,22,26,29,34,37,45),labels=1:6)

#gr_aci.a<-subset(gr_aci.a,gr_aci.a$TleafFac==2)

#--------------------------------------------------


tofit.2 <- split(gr_aci.a,gr_aci.a$Room)

g1fits<- get_topts(lapply(tofit.2,FUN=fitStom))
g1fits$Room<-c(1,4,6)
g1fits$Tgrowth<-c(18,28.5,35.5)



write.csv(g1fits,file = "Parameters/great_g1_ACi.csv",row.names=FALSE,sep=",")

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

#- estimate g1 from short term temperature response data

avt_short<-getAvT()
avt_short <- subset(avt_short, LightFac %in% c(4) & avt_short$Cond<2.5) #remove three outliers and get PAR>1500 

tofit.s <- split(avt_short,avt_short$Room)

g1fits<- get_topts(lapply(tofit.s,FUN=fitStom))
g1fits$Room<-c(seq(1,6,1))
g1fits$Tgrowth<-c(18,21.5,25,28.5,32,35.5)

write.csv(g1fits,file = "Parameters/great_g1_short_term.csv",row.names=FALSE,sep=",")

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------

# Interpolate leaf area, leaf mass, stem mass and root mass

# harvest_dat<-read.csv(paste0(path,"/Data/Glasshouse_DRAKE_EUTE_THERMAL-NICHE/data/harvest_data.csv"))


harvest_dat<-read.csv("Parameters/modelled_data.csv") # leaf area, leaf mass, stem mass and root mass from allometric models
harvest_dat$Date<-as.Date(harvest_dat$Date,format="%Y-%m-%d")

# harvest_dat_h<-read.csv("Parameters/modelled_data.csv") # get data from final harvest
# # harvest_dat_h<-subset(harvest_dat_h,harvest_dat_h$Date=="2016-02-22")
# harvest_dat_h$Date<-as.Date(harvest_dat_h$Date,format="%Y-%m-%d")

datesout <- seq(as.Date("2016/1/8"), as.Date("2016/2/29"), by="1 day") # get a sequance of dated 


# function to interpolate variables

interp_var <- function(d, datesout,var){
  
   sp <- spline(x=d$Date, y=d[[var]], xout=datesout) #smooth curve
   # sp <- approx(x=d$Date, y=d[[var]], xout=datesout) # linear interpolation
  
  out<-data.frame(Date=datesout, Room=unique(d$Room), var=sp$y)
  return(out)
}


#- leaf area
leaf_area <- do.call(rbind, lapply(split(harvest_dat, harvest_dat$Room), interp_var, datesout=datesout,var="Leafarea"))
names(leaf_area)[3]<-"leafarea_p"

#- leaf mass
Leafmass <- do.call(rbind, lapply(split(harvest_dat, harvest_dat$Room), interp_var, datesout=datesout,var="Leafmass"))
names(Leafmass)[3]<-"Leafmass"

#- stem mass

Stemmass <- do.call(rbind, lapply(split(harvest_dat, harvest_dat$Room), interp_var, datesout=datesout,var="Stemmass"))
names(Stemmass)[3]<-"Stemmass"

#- root mass
Rootmass <- do.call(rbind, lapply(split(harvest_dat, harvest_dat$Room), interp_var, datesout=datesout,var="Rootmass"))
names(Rootmass)[3]<-"Rootmass"


seedling_mass<-cbind(leaf_area,Leafmass[3],Stemmass[3],Rootmass[3])


write.csv(seedling_mass,file = "Parameters/great_leaf_area_predicted.csv",row.names=FALSE,sep=",")



#--Plot interpolated leaf area

COL <- rev(brewer.pal(6,"Spectral"))
windows(100,100);

par(cex.lab=1.5,mar=c(0,5,2,0.5),oma=c(2,0,0,0),cex.axis=1,las=1,mfrow=c(2,2))
plotBy(leafarea_p~Date|Room,type="l",col=COL,data=seedling_mass,legend=F,lwd=3,ylab="",xlab="")
with(harvest_dat,points(Date,Leafarea,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
# with(harvest_dat_h,points(Date,Leafarea,bg=COL[Room],pch=21,cex=2,ylab="",xlab="",axes=F))
# adderrorbars(x=harvest_dat_h$Date,y=harvest_dat_h$Leafarea,SE=harvest_dat_h$Leafarea_SE,direction="updown")

# axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)
#legend("topright",letters[1],bty="n",cex=1.2)

title(ylab=expression(Leaf[area]~(m^-2)),cex.lab=1.5,line=2.5)

legend("topleft",c("18","21.5","25","28.5","32.5","35.5"),fill=COL,cex=.8,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)


#--Plot leaf mass

plotBy(Leafmass~Date|Room,type="l",col=COL,data=seedling_mass,legend=F,lwd=3,ylab="",xlab="")
 with(harvest_dat,points(Date,Leafmass,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
title(ylab=expression(Leaf[mass]~(gC),cex.lab=1.5,line=1))

# with(harvest_dat_h,points(Date,Leafmass,bg=COL[Room],pch=21,cex=2,ylab="",xlab="",axes=F))
# adderrorbars(x=harvest_dat_h$Date,y=harvest_dat_h$Leafmass,SE=harvest_dat_h$Leafmass_SE,direction="updown")


#--Plot stem mass

plotBy(Stemmass~Date|Room,type="l",col=COL,data=seedling_mass,legend=F,lwd=3,ylab="",xlab="")
with(harvest_dat,points(Date,Stemmass,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
title(ylab=expression(Stem[mass]~(gC),cex.lab=1.5,line=2.5))
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)

# with(harvest_dat_h,points(Date,Stemmass,bg=COL[Room],pch=21,cex=2,ylab="",xlab="",axes=F))
# adderrorbars(x=harvest_dat_h$Date,y=harvest_dat_h$Stemmass,SE=harvest_dat_h$Stemmass_SE,direction="updown")

#--Plot Root mass

plotBy(Rootmass~Date|Room,type="l",col=COL,data=seedling_mass,legend=F,lwd=3,ylab="",xlab="")
with(harvest_dat,points(Date,Rootmass,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
title(ylab=expression(Root[mass]~(gC),cex.lab=1.5,line=3))
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)

# with(harvest_dat_h,points(Date,Rootmass,bg=COL[Room],pch=21,cex=2,ylab="",xlab="",axes=F))
# adderrorbars(x=harvest_dat_h$Date,y=harvest_dat_h$Rootmass,SE=harvest_dat_h$Rootmass_SE,direction="updown")

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#- This script process the leaf, stem, and root respiration data we measured at the end of the study
#  and makes two plots. (1) The actual respiration data measured at 25 deg C.
#  (2) Estiamted tissue-specific and whole-plant respiration rates at in-situ temperatures
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

#- get the data
Rdat <- returnRcomponents_new() 

# write.csv(Rdat,"/Data/Glasshouse_DRAKE_EUTE_THERMAL-NICHE/data/R_components.csv",row.names=FALSE,sep=",")

#- average across provenances in each room

Rdat.m <- summaryBy(Rarea+Rmass+Rmass_insitu+Rarea_insitu~Room+Tair+W_treatment+Organ,data=Rdat,FUN=c(mean,standard.error),na.rm=T)
#Rdat.m$Rarea.standard.error[which(is.na(Rdat.m$Rarea.standard.error))] <- 0




#- process data for in-situ estimates. Convert rates to mmol CO2 g-1 day-1
# Rdat$Rtotal_insitu <- NA
# Rdat$Rtotal_insitu[leaves] <- Rdat$Rmass_insitu[leaves]*Rdat$leafdm[leaves]/1000*60*60*24*1e-6  # all masses were recorded in mg
# Rdat$Rtotal_insitu[stems] <- Rdat$Rmass_insitu[stems]*Rdat$Stemmass[stems]/1000*60*60*24*1e-6
# Rdat$Rtotal_insitu[roots] <- Rdat$Rmass_insitu[roots]*Rdat$totalRoot[roots]/1000*60*60*24*1e-6

#- average across locations for well watered treatments.
Rdat_mean_insitu <- summaryBy(Rmass_insitu+Rarea_insitu~Tair+Organ,
                              FUN=c(mean,standard.error),data=subset(Rdat,W_treatment=="w"),na.rm=T)



#-- unit conversion nmol (CO2) g-1 (biomass) s-1 to g(C)g-1(biomass) day-1
Rdat_mean_insitu$Rmass_insitu.mean<-with(Rdat_mean_insitu,Rdat_mean_insitu$Rmass_insitu.mean*10^-9*60*60*24*12.0107/.48)
Rdat_mean_insitu$ Rmass_insitu.standard.error<-with(Rdat_mean_insitu,Rdat_mean_insitu$ Rmass_insitu.standard.error*10^-9*60*60*24*12.0107/.48)



#--reshape data

mn<-reshape2::dcast(Rdat_mean_insitu, Tair~ Organ,value.var=c("Rmass_insitu.mean"))
names(mn)[2:4]<-c("R_leaf","R_root","R_stem")

ses<-reshape2::dcast(Rdat_mean_insitu, Tair~ Organ,value.var=c("Rmass_insitu.standard.error"))
names(ses)[2:4]<-c("R_leaf_se","R_root_se","R_stem_se")

Rdat_mean_insitu_f<-merge(mn,ses,by="Tair")
Rdat_mean_insitu_f$Room<-seq(1,6,1)
Rdat_mean_insitu_f[1]<-NULL

# Rdat_mean_insitu_f$R_leaf<-ifelse(Rdat_mean_insitu_f$Room %in% c(5,6),Rdat_mean_insitu_f$R_leaf*.9,Rdat_mean_insitu_f$R_leaf)
# Rdat_mean_insitu_f$R_root<-ifelse(Rdat_mean_insitu_f$Room %in% c(5,6),Rdat_mean_insitu_f$R_root*.9,Rdat_mean_insitu_f$R_root)
# Rdat_mean_insitu_f$R_stem<-ifelse(Rdat_mean_insitu_f$Room %in% c(5,6),Rdat_mean_insitu_f$R_stem*.9,Rdat_mean_insitu_f$R_stem)


write.csv(Rdat_mean_insitu_f,file = "Parameters/great_Resp_leaf_shoot_root_fh.csv",row.names=FALSE,sep=",")


#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

# get leaf R at 25C (area basis to estimate leaf dark respiration; units: mumol m-2 s-1)

#- average leaf respiration (area basis) across locations for well watered treatments.
Rdat_mean_insitu <- summaryBy(Rarea+Rmass~Tair+Organ,
                              FUN=c(mean,standard.error),data=subset(Rdat,W_treatment=="w"),na.rm=T)

Rdat_mean_insitu<-subset(Rdat_mean_insitu,Rdat_mean_insitu$Organ=="leaf")[c(1,3:6)]


#-- unit conversion mumol (CO2) m-2 (LA) s-1 to gC m-2 s-1

# Rdat_mean_insitu$Rarea.mean<-with(Rdat_mean_insitu,Rdat_mean_insitu$Rarea.mean*10^-6*12.0107) # units gC m-2 s-1
# Rdat_mean_insitu$ Rarea.standard.error<-with(Rdat_mean_insitu,Rdat_mean_insitu$Rarea.standard.error*10^-6*12.0107)# units gC m-2 s-1

Rdat_mean_insitu$Room<-seq(1,6,1)
names(Rdat_mean_insitu)[2:5]<-c("R_leaf_a","R_leaf_m","R_leaf_a_se","R_leaf_m_se")
Rdat_mean_insitu$Tair<-NULL

write.csv(Rdat_mean_insitu,file = "Parameters/great_Resp_leaf_area_fh.csv",row.names=FALSE,sep=",")



# mn_d<-reshape2::dcast(Rdat_mean_insitu, Tair~ Organ,value.var=c("Rarea_insitu.mean"))
# names(mn_d)[2:4]<-c("R_leaf_a","R_root_a","R_stem_a")
# 
# ses_d<-reshape2::dcast(Rdat_mean_insitu, Tair~ Organ,value.var=c("Rarea_insitu.standard.error"))
# names(ses_d)[2:4]<-c("R_leaf_se","R_root_se","R_stem_se")
# 
# Rdat_mean_insitu_area<-merge(mn_d,ses_d,by="Tair")
# Rdat_mean_insitu_area$Room<-seq(1,6,1)
# Rdat_mean_insitu_area[c(1,3,4,6,7)]<-NULL
# 
# # Rdat_mean_insitu_area$Rd<-with(Rdat_mean_insitu_area,R_leaf_a*0.6)
# Rdat_mean_insitu_area$Rd<-with(Rdat_mean_insitu_area,R_leaf_a) # not corrected for Rdayfrac
# 
# write.csv(Rdat_mean_insitu_area,file = "Parameters/great_Resp_leaf_area_fh.csv",row.names=FALSE,sep=",")


#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------


Rdat_mean_insitu_f$Tgrowth<-c(18,21.5,25,28.5,32.5,35.5)
Rdat_mean_insitu_f$R_total<-with(Rdat_mean_insitu_f,R_leaf+R_root+R_stem)


#plot Rleaf
windows(100,100);
par(cex.lab=1.5,mar=c(4,7,2,0.5),oma=c(0,0,0,0),cex.axis=1.5,las=1,mfrow=c(2,2))
plot(R_leaf~Tgrowth,pch=1,col="red",data=Rdat_mean_insitu_f,legend=F,lwd=3,ylab="",xlab="",cex=2,ylim=c(0.02,0.07),axes=F)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
adderrorbars(x=Rdat_mean_insitu_f$Tgrowth,y=Rdat_mean_insitu_f$R_leaf,SE=Rdat_mean_insitu_f$R_leaf_se,direction="updown")
title(ylab=expression(R[leaf]~(gC~gC^-1~d^-1)),cex.lab=1.5,line=3.5)
title(xlab=expression(T[growth]~(degree*C),cex.lab=1.5,line=3))

# plot Root R
plot(R_root~Tgrowth,pch=1,col="red",data=Rdat_mean_insitu_f,legend=F,lwd=3,ylab="",xlab="",cex=2,ylim=c(0.02,0.07),axes=F)
adderrorbars(x=Rdat_mean_insitu_f$Tgrowth,y=Rdat_mean_insitu_f$R_root,SE=Rdat_mean_insitu_f$R_root_se,direction="updown")
title(ylab=expression(R[Root]~(gC~gC^-1~d^-1)),cex.lab=1.5,line=3.5)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
title(xlab=expression(T[growth]~(degree*C),cex.lab=1.5,line=3))

# plot Stem R
plot(R_stem~Tgrowth,pch=1,col="red",data=Rdat_mean_insitu_f,legend=F,lwd=3,ylab="",xlab="",cex=2,ylim=c(0.02,0.07),axes=F)
adderrorbars(x=Rdat_mean_insitu_f$Tgrowth,y=Rdat_mean_insitu_f$R_stem,SE=Rdat_mean_insitu_f$R_stem_se,direction="updown")
title(ylab=expression(R[Stem]~(gC~gC^-1~d^-1)),cex.lab=1.5,line=3.5)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)

title(xlab=expression(T[growth]~(degree*C),cex.lab=1.5,line=3))

# plot total R
plot(R_total~Tgrowth,pch=1,col="red",data=Rdat_mean_insitu_f,legend=F,lwd=3,ylab="",xlab="",cex=2,ylim=c(0.05,0.17),axes=F)
# adderrorbars(x=Rdat_mean_insitu_f$Tgrowth,y=Rdat_mean_insitu_f$R_stem,SE=Rdat_mean_insitu$R_stem_se,direction="updown")
title(ylab=expression(R[Total]~(gC~gC^-1~d^-1)),cex.lab=1.5,line=3.5)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)

title(xlab=expression(T[growth]~(degree*C),cex.lab=1.5,line=3))

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

#- get in-situ leaf respiration rates for from the short-term respiration data


Rdat <- data.frame(Room=1:6,Tair= c(18,21.5,25,28.5,32,35.5))

#- predict respiration at growth temperature
Rdat$Rmass_insitu <- with(Rdat,12.9*2^((Tair-22)/10)) #basal rate from Drake et al GCB 2017. mean of three provenances
Rdat$Rarea_insitu <- with(Rdat,0.68*2^((Tair-22)/10))
Rdat$Date<-as.Date("2016-02-11")

# unit conversion nmol (CO2) g-1 (biomass) s-1 to g(C)g-1(biomass) day-1
Rdat$Rmass_insitu<-with(Rdat,Rmass_insitu*10^-9*60*60*24*12.0107/.48)

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

# get stem and root respiration rates per second 

#- get the data
Rdat <- returnRcomponents() 

#--------------------------------------------------------------------------------------------------------------
#- average across locations for well watered treatments.
# Rdat_mean_insitu <- summaryBy(Rmass_insitu+Rarea_insitu~Tair+Organ,
#                               FUN=c(mean,standard.error),data=subset(Rdat,W_treatment=="w"),na.rm=T)

Rdat_mean_25 <- summaryBy(Rmass~Tair+Organ,
                          FUN=c(mean,standard.error),data=subset(Rdat,W_treatment=="w"),na.rm=T)

#-- unit conversion nmol (CO2) g-1 (biomass) s-1 to g(C)g-1(biomass) s-1

Rdat_mean_25$Rmass.mean<-with(Rdat_mean_25,Rmass.mean*10^-9*12.0107/.48)
Rdat_mean_25$Rmass.standard.error<-with(Rdat_mean_25,Rmass.standard.error*10^-9*12.0107/.48)



#--reshape data

mn<-reshape2::dcast(Rdat_mean_25, Tair~ Organ,value.var=c("Rmass.mean"))
names(mn)[2:4]<-c("R_leaf","R_root","R_stem")

ses<-reshape2::dcast(Rdat_mean_25, Tair~ Organ,value.var=c("Rmass.standard.error"))
names(ses)[2:4]<-c("R_leaf_se","R_root_se","R_stem_se")

Rdat_mean_insitu_f<-merge(mn,ses,by="Tair")
Rdat_mean_insitu_f$Room<-seq(1,6,1)
Rdat_mean_insitu_f[c(1)]<-NULL

# Rdat_mean_insitu_f$R_leaf<-ifelse(Rdat_mean_insitu_f$Room %in% c(5,6),Rdat_mean_insitu_f$R_leaf*.9,Rdat_mean_insitu_f$R_leaf)
# Rdat_mean_insitu_f$R_root<-ifelse(Rdat_mean_insitu_f$Room %in% c(5,6),Rdat_mean_insitu_f$R_root*.9,Rdat_mean_insitu_f$R_root)
# Rdat_mean_insitu_f$R_stem<-ifelse(Rdat_mean_insitu_f$Room %in% c(5,6),Rdat_mean_insitu_f$R_stem*.9,Rdat_mean_insitu_f$R_stem)


# get leaf R at 25C (area basis to estimate leaf dark respiration; units: mumol m-2 s-1)

#- average leaf respiration (area basis) across locations for well watered treatments.
Rdat_mean_leaf <- summaryBy(Rarea~Tair+Organ,
                              FUN=c(mean,standard.error),data=subset(Rdat,W_treatment=="w"),na.rm=T)

Rdat_mean_leaf<-subset(Rdat_mean_leaf,Rdat_mean_leaf$Organ=="leaf")
Rdat_mean_leaf[c(1,2)]<-NULL

Rdat_mean_leaf$Room<-seq(1,6,1)
names(Rdat_mean_leaf)[1:3]<-c("R_leaf_a","R_leaf_a_se","Room")


Rdat_mean_25<-merge(Rdat_mean_insitu_f,Rdat_mean_leaf,by="Room")

write.csv(Rdat_mean_25,file = "Parameters/great_Resp_leaf_shoot_root_25C.csv",row.names=FALSE,sep=",")







# Tgrowth<-c(18.0,22.5,25.0,28.5,32.5,35.5)
# Q10<-3.22-0.046*T_growth
# Q10_est<-data.frame(T_growth,Q10)






# #estimate alpha
# #alpha was estimated by linear regression at low PARi
# 
# #  Phooto=alpha[(PARi/4)*(Ci-Gammastar/Ci+Gammastar)]-RL
# 
# #Gammastar <- Bernachci et al 2001
# #Tresponse of Rday (RL): parameters from Drake et al 2017 GCB 
# 
# #- use the avt data (i.e., "long-term" photosynthesis dataset)
# avt_short<-getAQ()
# avt_short$GamStar<-TGammaStar(avt_short$Tleaf,Patm=avt_short$Press) #estimate Gammastar 
# 
# #get the independent variable for the model
# 
# avt_short$xt<-with(avt_short,PARi/4*((Ci-GamStar)/(Ci+2*GamStar)))
# 
# 
# #function to get the linear regression coefficients 
# 
# get_alpha<-function(data){
#   
#   lm_a<-lm(Photo~xt,data=data)
#   
#   alpha<-summary(lm_a)$coefficients[2]
#   alpha.se<-summary(lm_a)$coefficients[4]
#   param<-data.frame(cbind(alpha,alpha.se))
#   return(param)
#   
#   }
# avt_low<-subset(avt_short,avt_short$PARi<600 & avt_short$W_treatment=="w")
# avt_low_list<-split(avt_low,avt_low$Room)
# 
# alpha_est<-get_topts(lapply(avt_low_list,FUN=get_alpha))
# alpha_est$Room<-c(seq(1,6,1))
# alpha_est$Tgrowth<-c(18,21.5,25,28.5,32,35.5)
# 
# 
# write.csv(alpha_est,file = "Parameters/great_alpha_long_term.csv",row.names=FALSE,sep=",")
# #---------------------------------------------------------------------------------
# 
# 
# #----------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------
# 
# #estimate alpha with another method
# #alpha=4/PAR*(Ci+2*GammaStar)/(Ci-GammaStar)*(Photo+Rday)
# 
# #need some function for temperature response of Rday
# 
# #- get the avt data (i.e., "long-term" photosynthesis dataset)
# #avt <- getAvT()
# #avt$GamStar<-TGammaStar(avt$Tleaf,Patm=avt$Press) #estimate Gammastar 
# #avt$Rday<-TRday(avt$Tleaf)
# 
# avt_short$Rday<-TRdayQ10(Tleaf=avt_short$Tleaf)# Q10 model for day respiration; assumed as similer to Rdark
# 
# 
# avt_low<-subset(avt_short,avt_short$PARi<150 & avt_short$W_treatment=="w")
# avt_low$alpha<-with(avt_low,4/PARi*((Ci-2*GamStar)/(Ci-GamStar))*(Photo+Rday))
# 
# alpha_mean<-summaryBy(alpha~Room,data=avt_low,FUN=c(mean,std.error))
# alpha_mean$Tgrowth<-c(18,21.5,25,28.5,32,35.5)
# 
# write.csv(alpha_mean,file = "Parameters/great_alpha_long_term_invert_method.csv",row.names=FALSE,sep=",")
# 

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

# #estimate Vcmax, Jmax and their temperature response parameters
# 
#fit ACi curves
path<-getwd()
pathtodata<-paste0(path,"/Data/Glasshouse_DRAKE_EUTE_THERMAL-NICHE/data")


fg<-makecurves(pathtodata,"/GHS39_GREAT_MAIN_ACiT_20160216-20160227_L1.csv")
dg <- makedata.1(pathtodata,"/GHS39_GREAT_MAIN_ACiT_20160216-20160227_L1.csv",fg)
dg$TleafFac <-cut(dg$Tleaf,breaks=c(15,20,23,27,30,34,37),labels=1:6)
dg$JVr<-with(dg,Jmax/Vcmax)

 
#fit temperature response of Vcmax and Jmax
dg<-subset(dg,!dg$Curve %in% c(106,98)) #remove two outliers
dg.1<-split(dg,paste(dg$Room))
gr.vcmax<-get_topts(lapply(dg.1,function(x)fitvcmax_mm(x,random="Replicate",return="Peak")))

gr.vcmax$Room<-names(dg.1)
 
gr.vcmax[c(8:12)]<-NULL #remove some unwanted columns

 #fit temperature response of Jmax

 dg.a<-subset(dg,!is.na(Jmax) &!dg$Curve %in% c(82,96) ) #function do not work with NAs, so remove few NAs
 dg.2<-split(dg.a,paste(dg.a$Room))
 

 gr.jmax<-get_topts(lapply(dg.2,function(x)fitjmax_mm(x,random="Replicate",return="Peak")))
 
 gr.jmax$Room<-names(dg.2)

 gr.jmax[c(8:12)]<-NULL #remove some unwanted columns
 
 # to get vcmax and jmax in to one dataframe

 great_aci_params<-merge(gr.vcmax,gr.jmax,by="Room")
 great_aci_params$Tgrowth<-c(18,28.5,35.5)

 # #---------------------------------------------------------------------------------------------------------------
 # #---------------------------------------------------------------------------------------------------------------
 # 
 # #--- Vcmax and Jmax and their temperature resp[onse parameters are only available for 3 growth
 # #--- temperatures (18, 28.5 and 35.5). I used linear regression to estimate these parameters for 
 # #--- rest of the growth temperatures
 # 
 #--- function to fit smooth spline for vcmax and jmax
 
 interp_phy <- function(d, Tgrowth,var){
   
   sp <- spline(x=d$Tgrowth, y=d[[var]], xout=c(22.5,25,32.5))
   re <-data.frame(Tgrowth=Tgrowth, var=sp$y)
   return(re)
   
 }
 
 
 #--- get dataframe
 #--- only Vcmax25 and Jmax25 were estimated by smooth spline. Other parameters were set to 
 # #--- values at room 1 for rooms 2, room 4 for room 3 and values at room 6 for room 5
 # 
Tgrowth<-c(22.5,25,32.5)
Room<-c(2,3,5)
great_aci_rest<-data.frame(cbind(Tgrowth,Room))
 # 
 # #-- smooth curve
great_aci_rest$Jmax25<-interp_phy(great_aci_params,Tgrowth=c(22.5,25,32.5),var="Jmax25")[[2]]
great_aci_rest$Vcmax25<-interp_phy(great_aci_params,Tgrowth=c(22.5,25,32.5),var="Vcmax25")[[2]]

 
 
# assume temperature response not changed 

great_aci_rest$EaV<-great_aci_params$EaV
great_aci_rest$EaJ<-great_aci_params$EaJ
great_aci_rest$delsV<-great_aci_params$delsV
great_aci_rest$delsJ<-great_aci_params$delsJ

great_aci_all<-plyr::rbind.fill(great_aci_params,great_aci_rest)
write.csv(great_aci_all,file = "Parameters/temperature_response_parameters.csv",row.names=FALSE,sep=",")
 # 

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

# estimate Rday from Aci data

rday_dat<-summaryBy(Rd~Room,data=subset(dg,dg$TleafFac==3 & dg$Rd>0),FUN=c(mean,std.error))
rday_dat$Tgrowth<-c(18,25.8,35.5)
Tgrowth<-c(22.5,25,32.5)
Room<-c(2,3,5)
rday_dat_rest<-data.frame(cbind(Tgrowth,Room))

 
# #-- smooth curve
rday_dat_rest$Rd.mean<-interp_phy(rday_dat,Tgrowth=c(22.5,25,32.5),var="Rd.mean")[[2]]


great_rdayl<-plyr::rbind.fill(rday_dat_rest,rday_dat)
write.csv(great_rdayl[c(2,3)],file = "Parameters/rday_aci_fits.csv",row.names=FALSE,sep=",")

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------

# fit temperature response of Rday from ACi data

tofit<-subset(dg,dg$Rd>0)
tofit.a<-split(tofit,tofit$Room)
fitted_rday<-lapply(tofit.a,function(x)fitRvT(x,namex="Ts",namey="Rd"))

fitted_rday[[3]][[1]]
# JVr_mean<-summaryBy(JVr~TleafFac,data=dg,FUN=c(mean,std.error),na.rm=T)
# names(JVr_mean)[c(1:3)]<-c("Room",)
# 
# write.csv(JVr_sum,file = "Parameters/great_jvr_growth_t.csv",row.names=FALSE,sep=",")

# dg<-subset(dg,dg$Curve!=c(98)) #remove one outlier
# #dg.1<-split(dg,paste(dg$Room,dg$Provenance))
# 
# dg.1<-split(dg,paste(dg$Room))
# #dg.1<-split(dg,paste(dg$Provenance,dg$Room)) #no difference between provenances
# #dg.1[1]<-NULL #this is only four datapoints (ACi curves) whch has no identity (Room number or Provenance)
# 
# #dg.1<-split(dg,paste(dg$Provenance,dg$Room)) 
# #fit temperature response of Vcmax and Jmax
# 
# gr.vcmax<-get_topts(lapply(dg.1,function(x)fitvcmax_mm(x,random="Replicate",return="Peak")))
# #gr.vcmax$DataSet<-"GRATE"
# gr.vcmax$Room<-names(dg.1)
# #gr.vcmax$Provenance<-get_names(dg.1)[,1]
# gr.vcmax[c(8:12)]<-NULL #remove some unwanted columns
# 
# #fit temperature response of Jmax
# 
# dg.a<-subset(dg,!is.na(Jmax)) #function do not work with NAs, so remove few NAs
# dg.2<-split(dg.a,paste(dg.a$Room))
# #dg.2[1]<-NULL  #this is only four datapoints (ACi curves) whch has no identity (Room number or Provenance)
# 
# 
# gr.jmax<-get_topts(lapply(dg.2,function(x)fitjmax_mm(x,random="Replicate",return="Peak")))
# #gr.jmax$DataSet<-"GRATE"
# gr.jmax$Room<-names(dg.2)
# 
# gr.jmax[c(8:12)]<-NULL #remove some unwanted columns
# #gr.jmax$Provenance<-get_names(dg.2)[,1]
# 
# # to get vcmax and jmax in to one dataframe
# 
# great_aci_params<-merge(gr.vcmax,gr.jmax,by="Room")
# great_aci_params$Tgrowth<-c(18,28.5,35.5)
# 
# great_aci_params$Date<-as.Date("2016-02-16")
# 
# #----
# #----
# 
# # estimate J:V ratio at growth temperatures
# 
 JVr_sum<-summaryBy(JVr~Room+TleafFac,data=dg,fun=mean,keep.names=T,na.rm=T)
 # JVr_sum<-JVr_sum[c(1,9,17),]
 
 JVr_sum<-JVr_sum[c(3,10,17),]
 JVr_sum$Tgrowth<- c(18,28.5,35.5)
 write.csv(JVr_sum,file = "Parameters/great_jvr_25C.csv",row.names=FALSE,sep=",")
# 
#   
# # gr_vcmax_meas<-summaryBy(Vcmax+Jmax+Rd+Ts~Room+Provenance+TleafFac,data=subset(dg,!dg$Provenance %in% ""),FUN=mean,keep.names=T)
# # gr_vcmax_meas_b<-subset(gr_vcmax_meas, gr_vcmax_meas$Provenance=="B")
# with(dg,plot(Ts,Vcmax,col=Room,pch=16,cex=1.5))
# legend("topleft",legend=unique(dg$Room),col=unique(dg$Room),pch=19)
# 
# with(dg,plot(Ts,Jmax,col=Room,pch=16,cex=1.5))
# legend("topleft",legend=unique(dg$Room),col=unique(dg$Room),pch=19)
# 
# #---------------------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------------------

 #------------------------------------------------------------------------------------------------
 #------------------------------------------------------------------------------------------------
 
 #- get the data
 Rdat <- returnRcomponents_new() 
 
 #- average across provenances in each room
 
 Rdat.m <- summaryBy(Rarea+Rmass+Rmass_insitu+Rarea_insitu~Room+Tair+W_treatment+Organ,data=Rdat,FUN=c(mean,standard.error),na.rm=T)

 #- average across locations for well watered treatments.
 Rdat_mean_insitu <- summaryBy(Rmass_insitu+Rarea_insitu~Tair+Organ,
                               FUN=c(mean,standard.error),data=subset(Rdat,W_treatment=="w"),na.rm=T)
 
 
 
 #-- unit conversion nmol (CO2) g-1 (biomass) s-1 to g(C)g-1(biomass) day-1
 Rdat_mean_insitu$Rmass_insitu.mean<-with(Rdat_mean_insitu,Rdat_mean_insitu$Rmass_insitu.mean*10^-9*60*60*24*12.0107/.48)
 Rdat_mean_insitu$ Rmass_insitu.standard.error<-with(Rdat_mean_insitu,Rdat_mean_insitu$ Rmass_insitu.standard.error*10^-9*60*60*24*12.0107/.48)
 
 
 
 #--reshape data
 
 mn<-reshape2::dcast(Rdat_mean_insitu, Tair~ Organ,value.var=c("Rmass_insitu.mean"))
 names(mn)[2:4]<-c("R_leaf","R_root","R_stem")
 
 ses<-reshape2::dcast(Rdat_mean_insitu, Tair~ Organ,value.var=c("Rmass_insitu.standard.error"))
 names(ses)[2:4]<-c("R_leaf_se","R_root_se","R_stem_se")
 
 Rdat_mean_insitu_f<-merge(mn,ses,by="Tair")
 Rdat_mean_insitu_f$Room<-seq(1,6,1)
 Rdat_mean_insitu_f[1]<-NULL
 
 # Rdat_mean_insitu_f$R_leaf<-ifelse(Rdat_mean_insitu_f$Room %in% c(5,6),Rdat_mean_insitu_f$R_leaf*.9,Rdat_mean_insitu_f$R_leaf)
 # Rdat_mean_insitu_f$R_root<-ifelse(Rdat_mean_insitu_f$Room %in% c(5,6),Rdat_mean_insitu_f$R_root*.9,Rdat_mean_insitu_f$R_root)
 # Rdat_mean_insitu_f$R_stem<-ifelse(Rdat_mean_insitu_f$Room %in% c(5,6),Rdat_mean_insitu_f$R_stem*.9,Rdat_mean_insitu_f$R_stem)
 
 
 write.csv(Rdat_mean_insitu_f,file = "Parameters/great_Resp_leaf_shoot_root_fh.csv",row.names=FALSE,sep=",")
 
 
  
# #--- Vcmax and Jmax and their temperature resp[onse parameters are only available for 3 growth
# #--- temperatures (18, 28.5 and 35.5). I used linear regression to estimate these parameters for 
# #--- rest of the growth temperatures
# 
# #--- function to fit smooth spline for vcmax and jmax
#  
# interp_phy <- function(d, Tgrowth,var){
#   
#   sp <- spline(x=d$Tgrowth, y=d[[var]], xout=c(22.5,25,32.5))
#   re <-data.frame(Tgrowth=Tgrowth, var=sp$y)
#   return(re)
#   
# }
# 
# #--- get dataframe
# #--- only Vcmax25 and Jmax25 were estimated by smooth spline. Other parameters were set to 
# #--- values at room 1 for rooms 2, room 4 for room 3 and values at room 6 for room 5
# 
# Tgrowth<-c(22.5,25,32.5)
# Room<-c(2,3,5)
# 
# great_aci_rest<-data.frame(cbind(Tgrowth,Room))
# 
# #-- smooth curve
#   great_aci_rest$Jmax25<-interp_phy(great_aci_params,Tgrowth=c(22.5,25,32.5),var="Jmax25")[[2]]
#   great_aci_rest$Vcmax25<-interp_phy(great_aci_params,Tgrowth=c(22.5,25,32.5),var="Vcmax25")[[2]]
#   # 
# 
# # assume room 4 for room 2 & 3, room 6 for room 2
# # great_aci_rest$Vcmax25<-ifelse(great_aci_rest$Room==2,great_aci_params$Vcmax25[2],
# #                         ifelse(great_aci_rest$Room==3,great_aci_params$Vcmax25[2],great_aci_params$Vcmax25[3]))
# # 
# # great_aci_rest$Jmax25<-ifelse(great_aci_rest$Room==2,great_aci_params$Jmax25[2],
# #                                ifelse(great_aci_rest$Room==3,great_aci_params$Jmax25[2],great_aci_params$Jmax25[3]))
# 
# # assume room 1 for 2, room 4 for 3 and room 6 for 5
# great_aci_rest$Jmax25<-great_aci_params$Jmax25
# great_aci_rest$Vcmax25<-great_aci_params$Vcmax25
# 
# great_aci_rest$EaV<-great_aci_params$EaV
# great_aci_rest$EaJ<-great_aci_params$EaJ
# great_aci_rest$delsV<-great_aci_params$delsV
# great_aci_rest$delsJ<-great_aci_params$delsJ
# 
# great_aci_rest$Date<-as.Date("2016-02-16")
# great_aci_all<-plyr::rbind.fill(great_aci_params,great_aci_rest)
# 
# write.csv(great_aci_all,file = "yplantupscale/data/great_aci_t_response.csv",row.names=FALSE,sep=",") # for yplant
# write.csv(great_aci_all[c(1,3,4,10,11)],file = "Parameters/great_aci_t_response_fixed.csv",row.names=FALSE,sep=",") 
# 
# 
# #-------------------------------------------------------------
# 
# 
# 
# #-- get Vcmax and Jmax values at erlier time point
# #-- first, Vcmax and Jmax scalled to each room based on the photosynthesis ratio between two dates.
# #-- then, I assume Vcmax and Jmax of each room to be equal to the values of plants grown at 25C (room 3)
# #-- this values were assumed to be not changed from the start. so Vcmax and Jmax same for all rooms from 2016-01-07 to 2016-02-06
# #-- and decrease to measured values on 2016-02-16 and then no change untill harvest
# 
# #-- get Vcmax and Jmax values at erlier time point
# 
# table2 <- data.frame(Date=c(rep(as.Date("2016-02-03"),6)),Tgrowth=c(18,21.5,25,28.5,32.5,35.5),
#                      Room=seq(1:6))
# table2$Vcmax25<-NA
# table2$Vcmax25[which(table2$Room==1)]<-with(subset(great_aci_all,great_aci_all$Room==1),(Vcmax25-(Vcmax25*(0.009387171*(40-26)))))*.9
#  table2$Vcmax25[which(table2$Room ==2)]<-with(subset(great_aci_all,great_aci_all$Room ==2 ),(Vcmax25-(Vcmax25*(0.009387171*(40-26)))))*.9
# # table2$Vcmax25[which(table2$Room ==2)]<-with(subset(great_aci_all,great_aci_all$Room ==2 ),(Vcmax25-(Vcmax25*(-0.0302*(40-26)))))*.9
# 
# table2$Vcmax25[which(table2$Room ==3)]<-with(subset(great_aci_all,great_aci_all$Room ==3 ),(Vcmax25-(Vcmax25*(-0.0302*(40-26)))))*.9
# table2$Vcmax25[which(table2$Room ==4)]<-with(subset(great_aci_all,great_aci_all$Room ==4 ),(Vcmax25-(Vcmax25*(-0.0302*(40-26)))))*.9
# table2$Vcmax25[which(table2$Room ==5)]<-with(subset(great_aci_all,great_aci_all$Room ==5 ),(Vcmax25-(Vcmax25*(-0.0302*(40-26)))))*.9
# table2$Vcmax25[which(table2$Room ==6)]<-with(subset(great_aci_all,great_aci_all$Room ==6 ),(Vcmax25-(Vcmax25*(-0.0302*(40-26)))))*.9
# 
# 
# 
# table2$Jmax25<-NA
# table2$Jmax25[which(table2$Room==1)]<-with(subset(great_aci_all,great_aci_all$Room==1),(Jmax25-(Jmax25*(0.009387171*(40-26)))))*.9
#  table2$Jmax25[which(table2$Room ==2)]<-with(subset(great_aci_all,great_aci_all$Room ==2 ),(Jmax25-(Jmax25*(0.009387171*(40-26)))))*.9
# # table2$Jmax25[which(table2$Room ==2)]<-with(subset(great_aci_all,great_aci_all$Room ==2 ),(Jmax25-(Jmax25*(-0.0302*(40-26)))))*.9
# 
# table2$Jmax25[which(table2$Room ==3)]<-with(subset(great_aci_all,great_aci_all$Room ==3 ),(Jmax25-(Jmax25*(-0.0302*(40-26)))))*.9
# table2$Jmax25[which(table2$Room ==4)]<-with(subset(great_aci_all,great_aci_all$Room ==4 ),(Jmax25-(Jmax25*(-0.0302*(40-26)))))*.9
# table2$Jmax25[which(table2$Room ==5)]<-with(subset(great_aci_all,great_aci_all$Room ==5 ),(Jmax25-(Jmax25*(-0.0302*(40-26)))))*.9
# table2$Jmax25[which(table2$Room ==6)]<-with(subset(great_aci_all,great_aci_all$Room ==6 ),(Jmax25-(Jmax25*(-0.0302*(40-26)))))*.9
# 
# 
# 
# 
# #-- add T-response parameters (assume similar across dates)
# 
# table2<-merge(table2,great_aci_all[c(1,3,4,10,11)],by="Room")
# 
# 
# #-- get all data to a single table
# 
# great_aci_all<-rbind.fill(great_aci_all,table2)
# 
# 
# #-----------------------------------------------------------------------------
# #-- get the reduction in Vcmax and Jmax during 2016-02-03 to 2016-02-16
# #-- 
# 
# interp_var <- function(d, datesout,var){
#   
#   sp <- spline(x=d$Date, y=d[[var]], xout=datesout)
#   out<-data.frame(Date=datesout, Room=unique(d$Room), var=sp$y)
#   return(out)
# }
# 
# 
# datesout_1 <- seq(as.Date("2016/2/4"), as.Date("2016/2/15"), by="1 day")
# Vcmax_rest<-do.call(rbind, lapply(split(great_aci_all, great_aci_all$Room), interp_var, datesout=datesout_1,var="Vcmax25"))
# names(Vcmax_rest)[3]<-"Vcmax25"
# 
# Jmax_rest<-do.call(rbind, lapply(split(great_aci_all, great_aci_all$Room), interp_var, datesout=datesout_1,var="Jmax25"))
# names(Jmax_rest)[3]<-"Jmax25"
# 
# vc_rest<-merge(Vcmax_rest,Jmax_rest,by=c("Date","Room"))
# vc_rest<-merge(vc_rest,table2[c(1,3,6:9)],by="Room")
# 
# great_aci_params<-rbind.fill(great_aci_all,vc_rest)
# 



# get initial Vcmax and Jmax values. This was assumed as the Vcmax and Jmax values at Tgrowth=25C 
# on 2016-02-03 for all treatments

# Date<-as.Date("2016-01-07",format="%Y-%m-%d")
# Room<-c(1:6)
# Tgrowth<-c(18,22.5,25,28.5,32.5,35.5)
# table3<-data.frame(Date,Room,Tgrowth)
# 
# table3$Vcmax25<-great_aci_params[[2]][[9]]
# table3$Jmax25<-great_aci_params[[9]][[9]]
# 
# #- add T-response parameters
# table3<-merge(table3,table2[c(1,6:9)],by="Room")
# # table3<-merge(table3,table2[c(3,6,7,8)],by="Room")
# 
# great_aci_params<-rbind.fill(great_aci_params,table3)

write.csv(great_aci_params,file = "Parameters/great_aci_t_response.csv",row.names=FALSE,sep=",")

#----------------------------------------------------------------------------------------

# interp_phy <- function(d, Tgrowth,var){
#   
# sp <- spline(x=d$Tgrowth, y=d[[var]], xout=c(22.5,25,32.5))
# re <-data.frame(Tgrowth=Tgrowth, var=sp$y)
# return(re)
# 
# }
# 
# #sp <- spline(x=d$Tgrowth, y=d[["Vcmax25"]], xout=c(22.5,25,32.5))
# 
# 
# interp_phy(great_aci_params,Tgrowth=c(22.5,25,32.5),var="Jmax25")
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# #--- function to fit smooth spline for vcmax and jmax
# 
# interp_var_dates <- function(d,var){
#   
#   sp <- spline(x=d$Date, y=d[[var]], xout=seq(as.Date("2016/1/8"), as.Date("2016/3/1"), by = "day"))
#   re <-data.frame(Date=seq(as.Date("2016/1/8"), as.Date("2016/3/1"), by = "day"), var=sp$y)
#   return(re)
#   
# }



# #--- Get Final Drymass
# #- prepare biomass data, fit and plot the first panel
# dat <- getHarvest()
# size <- getSize()
# 
# #- pull out just the unique pot and room numbers from teh size dataframe
# size2 <- unique(size[,c("Code","Room","W_treatment","location","Tair")])
# 
# #- merge pot ids and harvest. Note the pre-treatment plants get excluded here
# dat2 <- merge(size2,dat,by.x=c("Code","location","W_treatment"),by.y=c("Code","location","W_treatment"))
# massdata <- subset(dat2,W_treatment == "w")
# 
# #- add TREATMENT MEANS for mass
# dat3 <- summaryBy(leafarea+leafdm+stemdm+rootdm+totdm+Tair~Room+Date,FUN=c(mean,standard.error),data=massdata,na.rm=T)
# dat3$Tair.standard.error<-NULL
# #-- add initial harvest
# 
# dat_in <- summaryBy(leafarea+leafdm+stemdm+rootdm+totdm~Room+Date,FUN=c(mean,standard.error),data=
#                     subset(dat,dat$Date==as.Date("2016-01-07")),na.rm=T)
# 
# dat_in<-dat_in[rep(seq_len(nrow(dat_in)), each=6),] # assume similar initial biomass for each room
# dat_in$Room<-c(1:6)
# dat_in$Tair.mean<-unique(dat3$Tair.mean)
# 
# 
# dat4<-rbind(dat_in,dat3)
# 
# #unit conversions (leaf area - m2)
# dat4$leafarea.mean<-dat4$leafarea.mean*10^-4
# dat4$leafarea.standard.error<-dat4$leafarea.standard.error*10^-4
# dat4$Tair.mean<-NULL
# names(dat4)<-c("Date","Leafarea","Leafmass","Stemmass","Rootmass","Totalmass","Leafarea_SE","Leafmass_SE",
#                "Stemmass_SE","Rootmass_SE","Totalmass_SE","Room")
# 
# write.csv(dat4,file = "Parameters/modelled_data.csv",row.names=FALSE,sep=",")


#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
windows()
COL.1=palette()[c(1:6)]
par(mfrow=c(2,2),mar=c(4,4,1,1),cex.axis=1.3,cex.lab=1.3)

with(dat4,plot(Date,Totalmass,col=COL.1[Room],pch=16,cex=2,ylab="Total mass (gC)"))
legend("topleft",legend=unique(dat4$Room),col=COL.1[dat4$Room],bty="n",pch=16,ncol=2,title="Room",cex=1.5)
adderrorbars(x=dat4$Date,y=dat4$Totalmass,SE=dat4$Totalmass_SE,direction="updown")
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)


with(dat4,plot(Date,Leafmass,col=COL.1[Room],pch=16,cex=2,ylab="Leaf mass (gC)"))
adderrorbars(x=dat4$Date,y=dat4$Leafmass,SE=dat4$Leafmass_SE,direction="updown")
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)

with(dat4,plot(Date,Stemmass,col=COL.1[Room],pch=16,cex=2,ylab="Stem mass (gC)"))
adderrorbars(x=dat4$Date,y=dat4$Stemmass,SE=dat4$Stemmass_SE,direction="updown")
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)


with(dat4,plot(Date,Rootmass,col=COL.1[Room],pch=16,cex=2,ylab="Root mass (gC)"))
adderrorbars(x=dat4$Date,y=dat4$Rootmass,SE=dat4$Rootmass_SE,direction="updown")
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)




# get harvest data 
#-----------------------------------------------------------------------------------------
################### Analyse stem height diameter to estimate Leaf, Stem and Root biomass
# Import weekly height diameter data for 3 months (Diameter is in mm; Height is in cm)
height.dia <- read.csv("Data/GHS39_GREAT_MAIN_HEIGHTDIAMETER_20160108-20160229_L2.csv")
height.dia = height.dia[height.dia$W_treatment %in% as.factor("w"),]
height.dia$D = rowMeans(height.dia[,c("D1", "D2")], na.rm=TRUE)
height.dia = height.dia[complete.cases(height.dia), ]
height.dia$Date = as.Date(height.dia$Date, format = "%d/%m/%Y")
# height.dia.sub = subset(height.dia, Date %in% as.factor("29/02/2016"))


initial.harvest = read.csv("Data/GHS39_GREAT_MAIN_BIOMASS_20160107_L2.csv")
initial.harvest[ , c("Leafmass", "Stemmass", "Rootmass")] = initial.harvest[ , c("Leafmass", "Stemmass", "Rootmass")]
initial.harvest$Height = initial.harvest$Height/10 # Unit conversion from mm to cm
initial.harvest$Date = as.Date("2016-01-07")
initial.harvest$Room = 0

int.harvest.1 = read.csv("Data/GHS39_GREAT_MAIN_BIOMASS_20160129_L2.csv")
int.harvest.1 = int.harvest.1[int.harvest.1$W_treatment %in% as.factor("w"),]
int.harvest.1[ , c("Leafmass", "Stemmass", "Rootmass")] = int.harvest.1[ , c("Leafmass", "Stemmass", "Rootmass")]
int.harvest.1$Date = as.Date("2016-01-29")
int.harvest.1 = unique(merge(int.harvest.1, height.dia[,c("Room","Pot")]))
int.harvest.2 = read.csv("Data/GHS39_GREAT_MAIN_BIOMASS_20160210_L2.csv")
int.harvest.2 = int.harvest.2[int.harvest.2$W_treatment %in% as.factor("w"),]
int.harvest.2[ , c("Leafmass", "Stemmass", "Rootmass")] = int.harvest.2[ , c("Leafmass", "Stemmass", "Rootmass")]
int.harvest.2$Date = as.Date("2016-02-10")
int.harvest.2 = unique(merge(int.harvest.2, height.dia[,c("Room","Pot")]))

final.harvest = read.csv("Data/GHS39_GREAT_MAIN_BIOMASS_20160217-20160224_L3_with date.csv")
final.harvest$Date = as.Date(final.harvest$Date, format = "%d/%m/%Y")
final.harvest = final.harvest[final.harvest$W_treatment %in% as.factor("w"),]
# final.harvest$Room = NULL
final.harvest$Leafarea = rowSums(final.harvest[,c("Leafarea", "Leafarea_sub")], na.rm=T)
final.harvest$Leafmass = rowSums(final.harvest[,c("Leafmass", "Leafmass_sub")], na.rm=T)
# final.harvest$Stemmass = rowSums(final.harvest[,c("Stemmass", "Stemmass_sub")], na.rm=T) # Clarify from John about the experiment procedure
final.harvest$Rootmass = rowSums(final.harvest[,c("Rootmass", "Rootmass_sub")], na.rm=T)
final.harvest[c("Leafarea_sub","Leafmass_sub","Stemmass_sub","Rootmass_sub")] = NULL
final.harvest[ , c("Leafmass", "Stemmass", "Rootmass")] = final.harvest[ , c("Leafmass", "Stemmass", "Rootmass")]/1000 # Unit conversion from mg to g
# final.harvest = unique(merge(final.harvest, height.dia[,c("Room","Pot")]))

harvest.data = rbind.fill(initial.harvest,int.harvest.1, int.harvest.2, final.harvest)
# harvest.data = rbind(int.harvest.1, int.harvest.2, subset(final.harvest, select = -Date))
harvest.data$W_treatment = NULL
# harvest.data = rbind(harvest.data, initial.harvest)
harvest.data$D = rowMeans(harvest.data[,c("D1", "D2")], na.rm=TRUE)
harvest.data = harvest.data[!(harvest.data$Rootmass==0),]
harvest.data = harvest.data[with(harvest.data, order(Rootmass)), ]

harvest.data$Date[which(harvest.data$Date %in% c(as.Date("2016-02-17"),as.Date("2016-02-18"),
                                                 as.Date("2016-02-19"),as.Date("2016-02-22"),as.Date("2016-02-23"),
                                                 as.Date("2016-02-24")))]<-as.Date("2016-02-22")

# harvest.data$Date[which(harvest.data$Date %in% c(as.Date("2016-02-22"),as.Date("2016-02-23"),
#                                                  as.Date("2016-02-24")))]<-as.Date("2016-02-24")



dat_mean<-summaryBy(Leafarea+Leafmass+Stemmass+Rootmass~Room+Date,data=harvest.data,FUN=c(mean,std.error))
names(dat_mean)<-c("Room","Date","Leafarea","Leafmass","Stemmass","Rootmass","Leafarea_SE","Leafmass_SE",
               "Stemmass_SE","Rootmass_SE")

dm_i<-rbind(dat_mean[1,],dat_mean[1,],dat_mean[1,],dat_mean[1,],dat_mean[1,],dat_mean[1,])
dm_i$Room<-seq(1,6,1)
dm_i$Date<-as.Date("2016-01-08")

data_harvest<-rbind(dat_mean,dm_i)
data_harvest<-subset(data_harvest,!data_harvest$Room==0)
data_harvest$Leafarea<-data_harvest$Leafarea*10^-4

data_harvest[ , c("Leafmass", "Stemmass", "Rootmass")] = data_harvest[ , c("Leafmass", "Stemmass", "Rootmass")]*0.48 # units in gC

write.csv(data_harvest,file = "Parameters/harvest_data.csv",row.names=FALSE,sep=",")

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
windows()
COL.1=palette()[c(1:6)]
par(mfrow=c(2,2),mar=c(4,4,1,1),cex.axis=1.3,cex.lab=1.3)

with(dat_mean,plot(Date,Totalmass,col=COL.1[Room],pch=16,cex=2,ylab="Total mass (gC)"))
legend("topleft",legend=unique(dat_mean$Room),col=COL.1[dat_mean$Room],bty="n",pch=16,ncol=2,title="Room")
adderrorbars(x=dat_mean$Date,y=dat_mean$Totalmass,SE=dat_mean$Totalmass_SE,direction="updown")
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)


with(dat_mean,plot(Date,Leafmass,col=COL.1[Room],pch=16,cex=2,ylab="Leaf mass (gC)"))
adderrorbars(x=dat_mean$Date,y=dat_mean$Leafmass,SE=dat_mean$Leafmass_SE,direction="updown")
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)

with(dat_mean,plot(Date,Stemmass,col=COL.1[Room],pch=16,cex=2,ylab="Stem mass (gC)"))
adderrorbars(x=dat_mean$Date,y=dat_mean$Stemmass,SE=dat_mean$Stemmass_SE,direction="updown")
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)


with(dat_mean,plot(Date,Rootmass,col=COL.1[Room],pch=16,cex=2,ylab="Root mass (gC)"))
adderrorbars(x=dat_mean$Date,y=dat_mean$Rootmass,SE=dat_mean$Rootmass_SE,direction="updown")
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)
legend("topleft",legend=unique(dat_mean$Room),col=COL.1[dat_mean$Room],bty="n",pch=16,ncol=2,title="Room")




#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------

# use original harvest data

# Interpolated leaf area, leaf mass, stem mass and root mass

harvest_dat<-read.csv(paste0(path,"/parameters/harvest_data.csv"))
harvest_dat$Date<-as.Date(harvest_dat$Date,format="%Y-%m-%d")

datesout <- seq(as.Date("2016/1/8"), as.Date("2016/2/29"), by="1 day")

interp_var <- function(d, datesout,var){
  
  # sp <- spline(x=d$Date, y=d[[var]], xout=datesout) #smooth curve
  sp <- approx(x=d$Date, y=d[[var]], xout=datesout) # linear interpollation
  
  out<-data.frame(Date=datesout, Room=unique(d$Room), var=sp$y)
  return(out)
}


leaf_area <- do.call(rbind, lapply(split(harvest_dat, harvest_dat$Room), interp_var, datesout=datesout,var="Leafarea"))
names(leaf_area)[3]<-"leafarea_p"

Leafmass <- do.call(rbind, lapply(split(harvest_dat, harvest_dat$Room), interp_var, datesout=datesout,var="Leafmass"))
names(Leafmass)[3]<-"Leafmass"

Stemmass <- do.call(rbind, lapply(split(harvest_dat, harvest_dat$Room), interp_var, datesout=datesout,var="Stemmass"))
names(Stemmass)[3]<-"Stemmass"

Rootmass <- do.call(rbind, lapply(split(harvest_dat, harvest_dat$Room), interp_var, datesout=datesout,var="Rootmass"))
names(Rootmass)[3]<-"Rootmass"


seedling_mass<-cbind(leaf_area,Leafmass[3],Stemmass[3],Rootmass[3])


write.csv(seedling_mass,file = "Parameters/great_leaf_area_predicted.csv",row.names=FALSE,sep=",") # used this for further calculations


#--Plot leaf area

COL <- rev(brewer.pal(6,"Spectral"))
windows(100,100);

par(cex.lab=1.5,mar=c(0,5,2,0.5),oma=c(2,0,0,0),cex.axis=1,las=1,mfrow=c(2,2))
plotBy(leafarea_p~Date|Room,type="l",col=COL,data=seedling_mass,legend=F,lwd=3,ylab="",xlab="")
with(harvest_dat,points(Date,Leafarea,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
# with(harvest_dat_h,points(Date,Leafarea,bg=COL[Room],pch=21,cex=2,ylab="",xlab="",axes=F))
# adderrorbars(x=harvest_dat$Date,y=harvest_dat$Leafarea,SE=harvest_dat$Leafarea_SE,direction="updown")

# axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)
#legend("topright",letters[1],bty="n",cex=1.2)

title(ylab=expression(Leaf[area]~(m^-2)),cex.lab=1.5,line=2.5)

legend("topleft",c("18","21.5","25","28.5","32.5","35.5"),fill=COL,cex=.8,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)


#--Plot leaf mass

plotBy(Leafmass~Date|Room,type="l",col=COL,data=seedling_mass,legend=F,lwd=3,ylab="",xlab="")
with(harvest_dat,points(Date,Leafmass,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
title(ylab=expression(Leaf[mass]~(gC),cex.lab=1.5,line=1))

# with(harvest_dat_h,points(Date,Leafmass,bg=COL[Room],pch=21,cex=2,ylab="",xlab="",axes=F))
# adderrorbars(x=harvest_dat_h$Date,y=harvest_dat_h$Leafmass,SE=harvest_dat_h$Leafmass_SE,direction="updown")


#--Plot stem mass

plotBy(Stemmass~Date|Room,type="l",col=COL,data=seedling_mass,legend=F,lwd=3,ylab="",xlab="")
with(harvest_dat,points(Date,Stemmass,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
title(ylab=expression(Stem[mass]~(gC),cex.lab=1.5,line=2.5))
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)

# with(harvest_dat_h,points(Date,Stemmass,bg=COL[Room],pch=21,cex=2,ylab="",xlab="",axes=F))
# adderrorbars(x=harvest_dat_h$Date,y=harvest_dat_h$Stemmass,SE=harvest_dat_h$Stemmass_SE,direction="updown")

#--Plot Root mass

plotBy(Rootmass~Date|Room,type="l",col=COL,data=seedling_mass,legend=F,lwd=3,ylab="",xlab="")
with(harvest_dat,points(Date,Rootmass,col=COL[Room],pch=16,cex=2,ylab="",xlab="",axes=F))
title(ylab=expression(Root[mass]~(gC),cex.lab=1.5,line=3))
axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-08"),to=as.Date("2016-03-10"),by="week"),labels=T)

# with(harvest_dat_h,points(Date,Rootmass,bg=COL[Room],pch=21,cex=2,ylab="",xlab="",axes=F))
# adderrorbars(x=harvest_dat_h$Date,y=harvest_dat_h$Rootmass,SE=harvest_dat_h$Rootmass_SE,direction="updown")










# with(subset(leaf_area,leaf_area$Room==1),plot(Date,leafarea_p,col=Room,pch=16,cex=1))
# with(subset(leaf_area,leaf_area$Room==6),points(Date,leafarea_p,col=Room,pch=16,cex=1))

# testR<-getRes(met,Rref_S=0.2,Rref_R=0.4,Q10=2.1)
# with(testR,plot(Tair,Rroot))
# testing purpose only: no need to run bellow code for main analysis
# #set up photosynthesis model model
# 
# #read GREAT climate data
# 
# met<-read.csv(paste0(pathtodata,"/GHS39_GREAT_MAIN_MET-AIR_20160107-20160302_L1.csv"))
# met<-subset(met,!is.na(PAR))
# 
# phy.data<-merge(gr.vcmax[c(1:3,14)],gr.jmax[c(1:3,14)],by="Room") #merge vcmax and jmax to one dataframe
# phy.data<-merge(g1fits[c(1,3)],phy.data,by="Room") #add g1
# phy.data<-merge(alpha_mean[c(1,2,4)],phy.data,by="Room") #add alpha
# 
# #metwithphy<-merge(met,phy.data,by="Room")
# 
# 
# #test simulation with AvT data
# 
# avtwithphy<-merge(avt,phy.data,by="Room")
# 
# #assume day respiration = respiration in light
# 
# photo_mod<-with(avtwithphy,Photosyn(VPD=VpdL, Ca=CO2S, PPFD=PARi, Tleaf=Tleaf,Patm=Press,
#                          gsmodel = "BBOpti", g1=g1,
#                          alpha = alpha.mean,
#                          theta = 0.85, Jmax =Jmax25, Vcmax =Vcmax25, Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
#                          EaV = EaV*10^3, EdVC = 2e+05, delsC = delsV*10^3,
#                          EaJ = EaJ*10^3, EdVJ = 2e+05, delsJ = delsJ*10^3))
# 
# with(avtwithphy,plot(Tleaf,Photo,col="green",pch=16,ylim=c(10,30)))
# with(photo_mod,points(Tleaf,ALEAF,col="red",pch=16))
# 
# 
# #################
# #room 1
# 
# room1<-subset(avtwithphy,avtwithphy$Room==1)
# photo_mod<-Photosyn(VPD=room1$VpdL, Ca=room1$CO2S, PPFD=room1$PARi, Tleaf=room1$Tleaf,Patm=room1$Press,
#                     gsmodel = "BBOpti",Jmax =210.65035, Vcmax =115.03459,EaV = 58.87154*10^3, EdVC = 2e+05, delsC = 0.6292704*10^3,
#                     EaJ = 42.74635*10^3, EdVJ = 2e+05, delsJ = 0.6312849*10^3,g1=1.578066,alpha=0.1927777)
# 
# with(room1,plot(Tleaf,Photo,col="green",pch=16,ylim=c(10,30)))
# with(photo_mod,points(Tleaf,ALEAF,col="red",pch=16))
# 
# 
# room4<-subset(avtwithphy,avtwithphy$Room==4)
# photo_mod_4<-Photosyn(VPD=room4$VpdL, Ca=room4$CO2S, PPFD=room4$PARi, Tleaf=room4$Tleaf,Patm=room4$Press,
#                     gsmodel = "BBOpti",Jmax =97.80397, Vcmax =69.06385,EaV = 57.71167*10^3, EdVC = 2e+05, delsC = 0.6288452*10^3,
#                     EaJ = 32.67277*10^3, EdVJ = 2e+05, delsJ = 0.6249740*10^3,g1=8.783680,alpha=0.2051992)
# 
# with(room4,plot(Tleaf,Photo,col="green",pch=16,ylim=c(10,40)))
# with(photo_mod_4,points(Tleaf,ALEAF,col="red",pch=16))
# 
# 
# ##################
# 
# #try to get photosynthesis numbers measured at ambient CO2 levels
# gr_aci<-read.csv(paste0(pathtodata,"/GHS39_GREAT_MAIN_ACiT_20160216-20160227_L0.csv"))
# gr_aci_amb<-subset(gr_aci,gr_aci$CO2R>430 & gr_aci$CO2R<450)
# 
# 
# aciwithphy<-merge(gr_aci_amb,phy.data,by="Room")
# #model photosynthesis for room 1, 4 and 6 using specified parameters for each room
# 
# #room 1
# photo_gr_1<-with(subset(aciwithphy,aciwithphy$Room==1),Photosyn(VPD=VpdL, Ca=CO2S, PPFD=PARi, Tleaf=Tleaf,Patm=Press,
#                                     gsmodel = "BBOpti", g1=g1,alpha = alpha.mean,EaV = EaV*10^3,EaJ = EaJ*10^3,
#                                     theta = 0.85, Jmax =Vcmax25*1.8, Vcmax =Vcmax25, Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
#                                     delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05
#                                     ))
# 
# with(subset(gr_aci_amb,gr_aci_amb$Room==1),plot(Tleaf,Photo,col="black",pch=16,ylim=c(0,40)))
# with(photo_gr_1,points(Tleaf,ALEAF,col="red",pch=16,ylim=c(0,40)))
# 
# 
# 
# #room 4
# 
# photo_gr_4<-with(subset(aciwithphy,aciwithphy$Room==4),Photosyn(VPD=VpdL, Ca=CO2S, PPFD=PARi, Tleaf=Tleaf,Patm=Press,
#                                                                 gsmodel = "BBOpti", g1=g1,
#                                                                 alpha = alpha.mean,
#                                                                 theta = 0.85, Jmax =Jmax25, Vcmax =Vcmax25, Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
#                                                                 EaV = EaV*10^3, EdVC = 2e+05, delsC = delsV*10^3,
#                                                                 EaJ = EaJ*10^3, EdVJ = 2e+05, delsJ = delsJ*10^3))
# 
# with(subset(gr_aci_amb,gr_aci_amb$Room==4),plot(Tleaf,Photo,col="black",pch=16,ylim=c(0,40)))
# with(photo_gr_4,points(Tleaf,ALEAF,col="red",pch=16,ylim=c(0,40)))
# 
# 
# 
# #room 6
# 
# photo_gr_6<-with(subset(aciwithphy,aciwithphy$Room==6),Photosyn(VPD=VpdL, Ca=CO2S, PPFD=PARi, Tleaf=Tleaf,Patm=Press,
#                                                                 gsmodel = "BBOpti", g1=g1,alpha = alpha.mean,
#                                                                 
#                                                                 theta = 0.85, Jmax =Jmax25, Vcmax =Vcmax25, Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
#                                                                 EaV = EaV*10^3, EdVC = 2e+05, delsC = delsV*10^3,
#                                                                 EaJ = EaJ*10^3, EdVJ = 2e+05, delsJ = delsJ*10^3))
# 
# with(subset(gr_aci_amb,gr_aci_amb$Room==6),plot(Tleaf,Photo,col="black",pch=16,ylim=c(0,40)))
# with(photo_gr_6,points(Tleaf,ALEAF,col="red",pch=16))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# tofit.l <- split(gr_aci_amb,paste(gr_aci_amb$Provenance,gr_aci_amb$Room))
# tofit.l[c(1)]<-NULL
# lapply(tofit.l,FUN=fitquad)
# 
# with(avtwithphy,plot(Tleaf,VpdL,col="red",pch=16))
# with(avtwithphy,plot(Tleaf,Cond,col="red",pch=16))
# with(avtwithphy,plot(Tleaf,Ci,col="red",pch=16))
# with(avtwithphy,plot(Tleaf,PARi,col="red",pch=16))
# with(avtwithphy,plot(Tleaf,Photo,col="red",pch=16))
# with(avtwithphy,plot(Tleaf,Trmmol,col="red",pch=16))
# with(avtwithphy,plot(Tleaf,CO2S,col="red",pch=16))
# 
# 
# avt$Photo_model<-with(avt,Cond/1.6*(CO2S-Ci))
# with(avt,plot(Photo,Photo_model,xlim=c(0,60),ylim=c(0,60)))
# abline(a=0,b=1)
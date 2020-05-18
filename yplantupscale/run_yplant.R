
#--- this script use YplantQMC R package to calculate self shading factor for different
#--- growth temperatures in GREAT experiment.

#--- 

# 
#  library(devtools)
#  install_bitbucket("remkoduursma/yplantqmc")

# if(!require(yplantqmc))install_bitbucket("remkoduursma/yplantqmc")
# pacman::p_load(YplantQMC, withr, lubridate, plantecophys)


if(!require(pacman))install.packages("pacman")
pacman::p_load(YplantQMC, withr, lubridate, plantecophys)

#devtools::install_bitbucket("remkoduursma/YplantQMC")
# Run this once (downloads quasimc executable)
#installQuasiMC()

# function to convert PPFD timeseries (for a day) to a total (MJ day-1)
par_to_total <- function(x)sum(x * 1E-06) * 15 * 60 / (4.57/2)

# Plants (3D information); these are all digitized Eucalyptus plants from Duursma et al. 2012 (NewPhyt)
euckey <- read.csv("yplantupscale/plantfiles/euc_key.csv", stringsAsFactors = FALSE)

eucs <- with_dir("yplantupscale/plantfiles", 
                 readplantlist(pfiles=euckey$pfile, lfiles=euckey$lfile))

# Calculate a summary for each plant, this calculates many variables.
eucsumm <- summary(eucs)

# Plot a plant (in 3D, can rotate plant)
#plot(eucs[[1]])


#--- set location with long and latitude (can plot() to see if it is correct)
richmond <- setLocation(lat=-33.6, long=150.75, tzlong=150)

#--- met data on 15 minute average from 2016-01-07 to 2016-03-02  (from GREAT met data)

weather <- read.csv("yplantupscale/data/s39_climate_15_min_averages.csv")
weather$DateTime_hr <- as.POSIXct(weather$DateTime_hr,format="%Y-%m-%d %T",tz="UTC")


weather <- transform(weather, 
                     DateTime15 = DateTime_hr,
                     Date = as.Date(DateTime_hr),
                     timeofday = hour(DateTime_hr) + minute(DateTime_hr)/60,
                     VPD = VPD)

# YplantQMC generates its own met data  
# Choose a day, calculate Tmin, Tmax, total PAR, maximum VPD, and use to generate met data.

#--- get a sunny day
sunday <- subset(weather, PAR >= 1800)
unique(sunday$Date)


sunny_date <- as.Date("2016-02-02") #same date as first in-situ gas exchange campaign
sunny_met <- subset(weather, Date == sunny_date)
sunny_met$DateTime_hr <- as.POSIXct(sunny_met$DateTime_hr,format="%Y-%m-%d %T",tz="UTC")

sunny_daymet <- setMet(richmond, month = month(sunny_date), day=day(sunny_date),
                 PARday = par_to_total(sunny_met$PAR), Tmin=min(sunny_met$Tair),
                 Tmax=max(sunny_met$Tair), VPDmax=max(sunny_met$VPD))

# to inspect,
plot(sunny_daymet)


#--- get a cloudy day

clday <- subset(weather,PAR<500)
unique(clday$Date)


cld_date <- as.Date("2016-01-27") #daily max PAR=390 mumols m-2s-1
cld_met <- subset(weather, Date == cld_date)
cld_daymet <- setMet(richmond, month = month(cld_date), day=day(cld_date),
                       PARday = par_to_total(cld_met$PAR), Tmin=min(cld_met$Tair),
                       Tmax=max(cld_met$Tair), VPDmax=max(cld_met$VPD))


# to inspect,
plot(cld_daymet)



# 3. Set physiological parameters for rooms seperately (set for room 1,4 and 6)

#--- Get Physiology parameters 

# grphy <- read.csv("yplantupscale/data/photo_physiology_parameters.csv")
grphy <- read.csv("yplantupscale/data/great_alpha_vcmax_jmax_g1.csv")

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
# read Rday (not corrected for rday frac =0.7)

rday<-read.csv("Parameters/great_Resp_leaf_area_fh.csv")

grphy<-merge(grphy,rday,by="Room")

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------





phy_to_run<-split(grphy,grphy$Room)  

# phy_to_run[c(2,3,5)]<-NULL #remove room 2,3 & 5


#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

#--- run yplant for a sunny day

room_list_sunny <- list()

for(i in 1:(length(phy_to_run))){
  
  d<-phy_to_run[[i]]
  # eucphy <- setPhy("Photosyn", leafpars=list(g1=d$g1,alpha = d$alpha.mean,EaV = d$EaV*10^3,EaJ = d$EaJ*10^3,
  #                                            theta = 0.85, Jmax =d$Jmax25, Vcmax =d$Vcmax25, Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
  #                                            delsC = d$delsV*10^3,delsJ = d$delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05))
  # 
  eucphy <- setPhy("Photosyn", leafpars=list(g1=d$g1,alpha = d$alpha,
                                             theta = 0.28, Jmax =d$Jmax25, Vcmax =d$Vcmax25,
                                             Tcorrect=F,Rd=d$Rd*0.6))
  
  
  room <- YplantDay(eucs, phy=eucphy, met=sunny_daymet)
  room.s<- summary(room)
  room.s$Room<-i
  room_list_sunny[[i]]<-room.s
  
}


saveRDS(room_list_sunny, "sunny.rds") 


yruns1 <- readRDS("sunny.rds")
sunny_summary<-do.call(rbind.data.frame,yruns1) 

#sunny_sum <- lapply(room_list_sunny, summary)

#need to add names the list by their room

#--------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------

#-- calculate self shading parameters for each growth temperature
#-- test linear regression to estimate self shading factor fron Leaf area

#-- add leaf area to simulations 
eucsumm_1<-eucsumm[,c("pfile", "LA")]
sunny_summary_la<-merge(sunny_summary,eucsumm_1,by="pfile")



#-- self shading factor
sunny_summary_la$self_s <- sunny_summary_la$totALEAF/sunny_summary_la$totALEAF0

# remove five out liers for testing
sunny_summary_la$self_s<-ifelse(sunny_summary_la$LA<1 & sunny_summary_la$self_s<0.88,0,sunny_summary_la$self_s)
sunny_summary_la<-subset(sunny_summary_la,sunny_summary_la$self_s>0)


mean_sigma<-summaryBy(self_s~Room,data=sunny_summary_la,FUN=c(mean,std.error))




#-- fit linear regression to get a relationship for self shading vs LA
library(broom)
sunny_list<-split(sunny_summary_la,sunny_summary_la$Room)
fit_sunny<-get_topts(lapply(sunny_list,function(x)fit_lm(x,"self_s","LA")))
fit_sunny$Room<-names(sunny_list)
# fit_sunny_param<-data.frame(do.call(rbind,list(fit_sunny[1][[1]],fit_sunny[2][[1]],fit_sunny[3][[1]],fit_sunny[4][[1]],
#                               fit_sunny[5][[1]],fit_sunny[6][[1]])))


fit_sunny<-merge(fit_sunny,mean_sigma,by="Room")


write.csv(fit_sunny,file = "Parameters/leaf_area_vs_self_shading.csv",row.names=FALSE,sep=",")

#--------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------

la_daily<-read.csv("Parameters/great_leaf_area_predicted.csv")
la_daily<-merge(la_daily,fit_sunny[c(1,2,4)],by="Room")
write.csv(la_daily,file = "Parameters/leaf_area_with_self_shading.csv",row.names=FALSE,sep=",")


# #--- run yplant for a cloudy day
# 
# room_list_cloudy <- list()
# 
# for(i in 1:(length(phy_to_run))){
#   
#   d<-phy_to_run[[i]]
#   eucphy <- setPhy("Photosyn", leafpars=list(g1=d$g1,alpha = d$alpha.mean,EaV = d$EaV*10^3,EaJ = d$EaJ*10^3,
#                                              theta = 0.85, Jmax =d$Jmax25, Vcmax =d$Vcmax25, Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
#                                              delsC = d$delsV*10^3,delsJ = d$delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05))
#   
#   room <- YplantDay(eucs, phy=eucphy, met=cld_daymet)
#   
#   room.s<- summary(room)
#   room.s$Room<-i
#   room_list_cloudy[[i]]<-room.s
#   
#   
#   
#   room_list_cloudy[[i]]<-room
#   
# }
# 
# saveRDS(room_list_sunny, "cloudy.rds")
# cloudy_summary<-do.call(rbind.data.frame,room_list_cloudy)

#--------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------

# room_test <- YplantDay(eucs[[1]], phy=eucphy, met=sunny_daymet)
# psrdata(room_test)
# plot(room_test)
# # Plot simulation results for one of the plants:
# # plot(yruns1[[1]])
# # plot(yruns2[[1]])
# 
# yres1 <- cbind(sunny_summary, eucsumm)
# yres2 <- cbind(cloudy_summary, eucsumm)
# 
# 
# with(yres1, plot(nleavesp, totA / totA0, ylim=c(0.7,1)))
# with(yres2, points(nleavesp, totA / totA0, pch=19))
# 


# sunny_summary_la$self_s<-ifelse(sunny_summary_la$LA<1 & sunny_summary_la$self_s<0.88,0,sunny_summary_la$self_s)
# sunny_summary_la<-subset(sunny_summary_la,sunny_summary_la$self_s>0)
# 
# with(subset(sunny_summary_la,sunny_summary_la$Room==1),plot(LA,self_s,col=Room,pch=16))
# with(subset(sunny_summary_la,sunny_summary_la$Room==2),plot(LA,self_s,col=Room,pch=16))
# with(subset(sunny_summary_la,sunny_summary_la$Room==3),plot(LA,self_s,col=Room,pch=16))
# with(subset(sunny_summary_la,sunny_summary_la$Room==4),plot(LA,self_s,col=Room,pch=16))
# with(subset(sunny_summary_la,sunny_summary_la$Room==5),plot(LA,self_s,col=Room,pch=16))
# with(subset(sunny_summary_la,sunny_summary_la$Room==6),plot(LA,self_s,col=Room,pch=16))


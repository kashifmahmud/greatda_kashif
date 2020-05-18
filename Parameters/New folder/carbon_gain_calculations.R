#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

#--- get required functions

# source("R/functions_for_analysis.R")
# source("R/generic_functions.R")

library(plantecophys)
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

#--- Daily carbon gain calculations for GREAT experiments


#--- read met data and fitted parameters

#--- read GREAT climate data

met<-read.csv("Parameters/s39_climate_15_min_averages.csv")

#--- read vcmax, jmax and temperature response parameters

bioch<-read.csv("Parameters/great_aci_t_response.csv")


#--- read g1

short_g1<-read.csv("Parameters/great_g1_short_term.csv") # short-term tempearture response


#--- read alpha

alph<-read.csv("Parameters/great_alpha_long_term_invert_method.csv") # alpha from invert method


#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# make dataframe with g1 and alpha

phy.data<-merge(short_g1[c(1:3)],alph[c(1:3)],by="Room") 

#add biochemical parameters 

phy.data<- merge(phy.data,bioch[c(1:4,9:11,16)],by="Room",all=T) 


#--- add physiological parameters to med data

metwithphy<-merge(met,phy.data,by="Room")

#apply(room3, 2, function(x) any(is.na(x)))

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

#set up photosynthesis model 

#Day respiration: Basal rate (at 22C) and Q10 were assumed to be simila to dark respiration rate
#Values from Drake et al 2017 GCB; averaged across provenances

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------



#make negatine PAR values (recorded in night) to zero.

metwithphy$PAR<-ifelse(metwithphy$PAR<0,0,metwithphy$PAR)
metwithphy<-subset(metwithphy,!is.na(PAR))

# for all rooms together (this is not run due to the error message)
photo_gr_all<-with(metwithphy,Photosyn(VPD=VPD, Ca=400,Tleaf=Tair,PPFD=PAR,gsmodel = "BBOpti", g1=g1,alpha = alpha.mean,EaV = EaV*10^3,EaJ = EaJ*10^3,
                                                                         theta = 0.85, Jmax =Jmax25, Vcmax =Vcmax25, Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
                                                                         delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05))


# room 1,2,5,6 together (no error)

photo_gr<-with(subset(metwithphy,metwithphy$Room %in% c(1,2,5,6)),Photosyn(VPD=VPD, Ca=400, Tleaf=Tair,PPFD=PAR, 
                                                                gsmodel = "BBOpti", g1=g1,alpha = alpha.mean,EaV = EaV*10^3,EaJ = EaJ*10^3,
                                                                theta = 0.85, Jmax =Jmax25, Vcmax =Vcmax25, Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
                                                                delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05))

photom.1 <- photo_gr[,c(1:4, 7:11)]
photom_15min_1 <- cbind(photom.1,subset(metwithphy,metwithphy$Room %in% c(1,2,5,6)))



# room 3,4 together (no error when Ca set to 399 ppm)
photo_gr.1<-with(subset(metwithphy,metwithphy$Room %in% c(3,4)),Photosyn(VPD=VPD, Ca=399, Tleaf=Tair,PPFD=PAR, 
                                                                           gsmodel = "BBOpti", g1=g1,alpha = alpha.mean,EaV = EaV*10^3,EaJ = EaJ*10^3,
                                                                           theta = 0.85, Jmax =Jmax25, Vcmax =Vcmax25, Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
                                                                           delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05))

photom.2 <- photo_gr.1[,c(1:4, 7:11)]
photom_15min_2 <- cbind(photom.2,subset(metwithphy,metwithphy$Room %in% c(3,4)))


photom_15min<-rbind(photom_15min_1,photom_15min_2)


# photom <- photo_gr[,c(1:4, 7:11)]
# photom_15min <- cbind(photom,subset(metwithphy,metwithphy$Room %in% c(1,2,5,6)))

photom_15min$Date <- as.Date(photom_15min$Date)

photom_15min$gCarbon <- with(photom_15min, ALEAF*15*60*10^-6*12)

dCarbon <- summaryBy(gCarbon ~ Date+Room, data=photom_15min, FUN=sum, keep.names=TRUE )
#names(dCarbon)[3] <- "carbon_day"
with(dCarbon,plot(Date,gCarbon,col=Room,pch=16))



write.csv(Aleaf_15min, "calculated data/daily_carbon_gain.csv", row.names=FALSE)

#--------------

photo_gr.1<-with(subset(metwithphy,metwithphy$Room %in% c(3,4)),Photosyn(VPD=VPD, Ca=399, Tleaf=Tair,PPFD=PAR, 
                                                                         gsmodel = "BBOpti", g1=g1,alpha = alpha.mean,EaV = EaV*10^3,EaJ = EaJ*10^3,
                                                                         theta = 0.85, Jmax =Jmax25, Vcmax =Vcmax25, Rd0 = 0.68, Q10 = 2.1,TrefR = 22,
                                                                         delsC = delsV*10^3,delsJ = delsJ*10^3,EdVC = 2e+05,EdVJ = 2e+05))


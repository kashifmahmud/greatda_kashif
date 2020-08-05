#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------


#--- read met data and fitted parameters

#--- read GREAT met data 15 min averages

met<-read.csv("Data/MET_DATA/s39_climate_15_min_averages.csv")
met$Date<-as.Date(met$Date)

gpp_na<-subset(met,met$Date %in% c(as.Date("2016-1-18"),as.Date("2016-1-19")))
gpp_20<-subset(met,met$Date %in% as.Date("2016-1-19"))
gpp_20$Date<-as.Date("2016-01-20")
gpp_21<-subset(met,met$Date %in% as.Date("2016-1-22"))
gpp_21$Date<-as.Date("2016-01-21")

# gpp_na$Date<-ifelse(gpp_na$Date=="2016-1-18","2016-1-20","2016-1-21")  
met<-rbind(met,rbind(gpp_20,gpp_21)) 

# met$Tair<-ifelse(met$Room %in% c(4,5),met$Tair-2,met$Tair)

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------



# Calculate respiration rates

resp_dat<-read.csv("parameters/great_Resp_leaf_shoot_root_25C.csv")
met<-merge(met,resp_dat,by=c("Room"))

#unit conversion of respiration rate at 25C (gC gC-1 s-1 to gC gC-1 15min)

met$R_stem<-with(met,R_stem*15*60)
met$R_root<-with(met,R_root*15*60)
met$R_leaf<-with(met,R_leaf*15*60)


met$Rdark_stem<-with(met, R_stem*2.1^((Tair-25)/10)) # units: gC gC-1 15min-1
met$Rdark_root<- with(met, R_root*2.1^((Tair-25)/10))# units: gC gC-1 15min-1
met$Rdark_leaf<- with(met, R_leaf*2.1^((Tair-25)/10))# units: gC gC-1 15min-1

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------


#calculate daily values (sum across dates)

met_sum<-summaryBy(Rdark_stem+Rdark_root+Rdark_leaf~Room+Date,data=met,FUN=sum,keep.names=T) #units gC gC-1 d-1

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
# Read tree attributes data including Met data, Temperature dependant variables, Modelled Parameters from Dushan's analysis 
data.attrib <- read.csv("parameters/data_for_attribution_analysis_v2.csv")
data.attrib$Date = as.Date(data.attrib$Date, format="%d/%m/%Y")
data.attrib$DateTime_hr = as.POSIXct(as.character(data.attrib$DateTime_hr), format="%d/%m/%Y %H:%M")

data.attrib.rm1 = subset(data.attrib, Room %in% c(1))
data.attrib.rm1 = subset(data.attrib.rm1, Date >= as.Date("2016-01-17") & Date <= as.Date("2016-01-22"))
par(mfrow=c(2,1))
plot(data.attrib.rm1$DateTime_hr, data.attrib.rm1$Tair, type = 'l')

met$DateTime_hr = as.POSIXct(as.character(met$DateTime_hr), format="%Y-%m-%d %H:%M:%S")
met.rm1 = subset(met, Room %in% c(1))
met.rm1 = subset(met.rm1, Date >= as.Date("2016-01-17") & Date <= as.Date("2016-01-22"))
met.rm1 = met.rm1[ order(met.rm1$DateTime_hr , decreasing = FALSE ),]
plot(met.rm1$DateTime_hr, met.rm1$Tair, type = 'l')



#---Read and process Aspinwall et al 2016 ACi data


path<-getwd()
#--- Read data

mike<-read.csv(paste0(path,"/Data/WTC_TEMP_CM_GX-ACI25_20130705-20140422_L1_v1.csv"))
mike<-subset(mike,mike$Tofit==1)

#add curev number
mike$Curve <- c(1)
Curve <- c()
count <- 1  
for (i in 2:length(mike$Code)){
  
  ifelse(mike$Code[i-1] == mike$Code[i],count <- count,count <- count + 1)
  Curve[i] <- count 
}

mike$Curve[2:length(mike$Code)] <- na.omit(Curve)
mike$Date <- as.POSIXct(mike$Date,format="%Y-%m-%d")

mike$age_of_plant<-as.numeric(with(mike,(Date-as.POSIXlt("2013-03-12")))) #no of days after planted inside WTC


#-- change treatment names (capital first letter)

capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}

mike$T_treatment<-capFirst(mike$T_treatment)
write.csv(mike,file="Data/WTC_TEMP_CM_GX-ACI25_20130705-20140422_L2.csv",row.names=F,sep=",")


#---Fit ACi curves

#fit ACi curves

pathtodata<-paste0(path,"/Data")

ma<-makecurves(pathtodata,"/WTC_TEMP_CM_GX-ACI25_20130705-20140422_L2.csv")
dma <- makedata.1(pathtodata,"/WTC_TEMP_CM_GX-ACI25_20130705-20140422_L2.csv",ma)
dma$Date <- as.POSIXct(dma$Date,format="%Y-%m-%d")

#----
#----


#--- function to get growth temperatures for ambient and warmed treatments (30 days presceding temperature)

wtcMet<- function(path,fname,from){
  
  #set WD to the data directory
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  #read all csv files and bind to one file
  
  wtcmet <- read.csv(paste0(path,"/",fname),stringsAsFactors=FALSE)
   
  wtcmet$DateTime <- as.POSIXct(wtcmet$DateTime,format="%Y-%m-%d")
  
  from=as.POSIXct(from,format="%Y-%m-%d")
  
  end30<-from - as.difftime(30, unit="days")
  # end60<-as.Date(from-15)
  # end90<-as.Date(from-3)
  # 
  #get a subset for required period
  dat30<-wtcmet[wtcmet$DateTime >= end30 & wtcmet$DateTime <= from,]
  # dat60<-wtcmet[wtcmet$DateTime >= end60 & wtcmet$DateTime <= from,]
  # dat90<-wtcmet[wtcmet$DateTime >= end90 & wtcmet$DateTime <= from,]
  # 
  
  
  #get average across given time period (set only chamber temperature)
  
  t30 <- dplyr::summarize(group_by(dat30,T_treatment),
                          #Tair_SP=mean(Tair_SP,na.rm=T),
                          #RH_al=mean(RH_al,na.rm=T),
                          #DP_al=mean(DP_al,na.rm=T),
                          Tmin=mean(Tmin,na.rm=T),
                          Tmax=mean(Tmax,na.rm=T),
                          Tmean=mean(Tmean,na.rm=T))
  
  
  # t60 <- dplyr::summarize(group_by(dat60,Treat),
  #                         #Tair_SP=mean(Tair_SP,na.rm=T),
  #                         #RH_al=mean(RH_al,na.rm=T),
  #                         #DP_al=mean(DP_al,na.rm=T),
  #                         Tmin=mean(Tmin,na.rm=T),
  #                         Tmax=mean(Tmax,na.rm=T),
  #                         Tmean=mean(Tmean,na.rm=T))
  # 
  # t90 <- dplyr::summarize(group_by(dat90,Treat),
  #                         #Tair_SP=mean(Tair_SP,na.rm=T),
  #                         #RH_al=mean(RH_al,na.rm=T),
  #                         #DP_al=mean(DP_al,na.rm=T),
  #                         Tmin=mean(Tmin,na.rm=T),
  #                         Tmax=mean(Tmax,na.rm=T),
  #                         Tmean=mean(Tmean,na.rm=T))
  # 
  
  taverage <- data.frame(t30)
  taverage$Date<-as.POSIXct(from,format="%Y-%m-%d")
  
  
  return(taverage)
}


#-- growth temperatures for measurement dates

met_dat<-list()
dates<-unique(dma$Date)
for(i in 1:length(dates)){
  
  from<-dates[i]
  met<-wtcMet(path=pathtodata,fname="WTC_TEMP_CM_WTCMET_L1_v2.csv",from=from)
  met_dat[[i]]<-met
}

met_data<-do.call(rbind.data.frame,met_dat)
names(met_data)[1]<-"T_treatment"

#--- add met data to the Vcmax dataframe

dma_with_met<-merge(dma,met_data,by=c("Date","T_treatment"))
dma_with_met<-subset(dma_with_met,!dma_with_met$Curve %in% c(37,38,71)) # remove 3 outlier datapoints


#---- plot Vcmax and Jmax vs age of seedlings

wtc_mean<-summaryBy(Vcmax+Jmax+age_of_plant+Tmean~Date+T_treatment,FUN=c(mean,std.error),data=dma_with_met)

with(wtc_mean,plot(age_of_plant.mean,Vcmax.mean,col=T_treatment,pch=16,cex=2))
with(wtc_mean,plot(Tmean.mean,Vcmax.mean,col=T_treatment,pch=16,cex=2))



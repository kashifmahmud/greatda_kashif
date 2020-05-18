

# library(plantecophys)
# 
# # get data at low PAR and split by room
# dat <- read.csv("C:/GreatGrowth_NA/Data/Glasshouse_DRAKE_EUTE_THERMAL-NICHE/data/GHS39_GREAT_MAIN_GX-AQ_20160202-20160203_L3.csv")

avt_long<-getAQ()
dat <- subset(avt_long,avt_long$W_treatment=="w" & avt_long$PARi<150)

# dat <- subset(dat,PARi < 150)
# dat$Room[which(dat$Room==4)]<-2
# dat$Room[which(dat$Room==6)]<-3

#-- add physiological data to the dataframe


#--- read g1

# short_g1<-read.csv("Parameters/great_g1_long_term.csv") # short-term tempearture response

short_g1<-read.csv("Parameters/great_g1_short_term.csv") # short-term tempearture response


#--- read ACi parameter

# aci_param<-read.csv("yplantupscale/data/great_aci_t_response.csv")[c(1:4,9:11)] 
aci_param<-read.csv("Parameters/temperature_response_parameters.csv")[c(1:4,9:11)]


#---

dat_ph<-merge(short_g1,aci_param,by="Room")



#---

dat_all<-merge(dat,dat_ph,by="Room")


datlist <- split(dat_all,dat_all$Room)

# function calling Photosyn
# DUSHAN - you will need to pass here the other parameters that you need to pass to Photosyn
# e.g. g1, delSj, Eaj, Rd, Q10 etc
pslow <- function(alpha,PAR,Tleaf,Vcmax25,Jmax25,EaV,EaJ,delsV,delsJ,EdVC=NULL,Rd0=NULL,q10=NULL,TrefR=NULL,Rdayfrac=NULL,g1,gsmodel=NULL,
                  VPD,CO2S) {
  
  ps <- Photosyn(alpha=alpha,PPFD=PAR,Tleaf=Tleaf,Vcmax=Vcmax25,Jmax=Jmax25,EaV=EaV*10^3,EaJ=EaJ*10^3,delsC=delsV*10^3,delsJ=delsJ*10^3,EdVC=2e+05,
                 Rd0=0.68,Q10=2.1,TrefR=22,Rdayfrac=0.7,g1=g1,gsmodel="BBOpti",VPD=VPD,Ca=CO2S)
  ret <- ps$Aj - ps$Rd
  return(ret)
}

# check function
# pslow(0.2,100,25,50,100)

# call to NLS
alpha<-alpha.se <- theta<-theta.se<-room <- tair <- c()

for (i in 1:length(unique(dat$Room))) {
  fitalpha <- nls(Photo ~ pslow(alpha,PAR=PARi, Tleaf=Tleaf,Vcmax=Vcmax25,Jmax=Jmax25,EaV=EaV,EaJ=EaJ,delsV=delsV,delsJ=delsJ,g1=g1,VPD=VpdL,CO2S=CO2S),
             start=list(alpha=0.3), data=datlist[[i]], trace=TRUE)
  
  alpha[i]<-summary(fitalpha)$parameters[[1]]
  # theta[i]<-summary(fitalpha)$parameters[[2]]
  alpha.se[i]<-summary(fitalpha)$parameters[[2]]
  # theta.se[i]<-summary(fitalpha)$parameters[[4]]
  room[i] <- i
  tair[i] <- mean(datlist[[i]]$Tair)
}

# look at output - ordered by room
alpha_est<-data.frame(room,tair,alpha,alpha.se)
names(alpha_est)[1]<-"Room"
write.csv(alpha_est[c(1,3,4)],file = "Parameters/great_alpha_long_term_nls_method.csv",row.names=FALSE,sep=",")



# plots
# modelled vs measured
plot(dat_all$Photo,pslow(alpha[dat_all$Room],dat_all$PARi,dat_all$Tleaf,dat_all$Vcmax25,dat_all$Jmax25,dat_all$EaV,dat_all$EaJ,dat_all$delsV,dat_all$delsJ,g1=dat_all$g1,VPD=dat_all$VpdL,CO2S=dat_all$CO2S),col=dat_all$Room,xlab="",ylab="")
abline(0,1)
# both vs temperature
with(dat_all,plot(Tleaf,Photo,col=Room))
with(dat_all,points(Tleaf,pslow(alpha[dat_all$Room],dat_all$PARi,dat_all$Tleaf,dat_all$Vcmax25,dat_all$Jmax25,dat_all$EaV,dat_all$EaJ,dat_all$delsV,dat_all$delsJ,g1=dat_all$g1,VPD=dat_all$VpdL,CO2S=dat_all$CO2S),col=dat_all$Room,pch=19))
legend("topright",legend=unique(dat_all$Room),col=unique(dat_all$Room),pch=19)

#---plot alpha vs Tgrowth

with(alpha_est,plot(tair,alpha,ylim=c(0.2,0.4)))




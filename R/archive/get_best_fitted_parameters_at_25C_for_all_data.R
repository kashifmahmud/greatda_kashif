
# code to get best fitted Vcmax and alpha for GREAT-AQ data
# temperature response parameters of Vcmax and Jmax: from ACi data
# day respiration: assume similar to Rdark (short term response curves) with Rdayfrac=0.7
# g1 estimated from AQ data Campaign 1 (In-situ photosynthesis measurements)

#------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------


# function calling Photosyn
# return Vcmax (at growth temperature) and alpha
# Jmax = J:V*Vcmax 
# J:V ratio estimated from ACi data

# function to call Photosyn

# pslow <- function(alpha,PAR,theta,Tleaf,Vcmax,Jmax,Rd0=NULL,q10=NULL,TrefR=NULL,Rdayfrac=NULL,g1,gsmodel=NULL,
#                   VPD,CO2S) {
# 
#   ps <- Photosyn(alpha=alpha,theta=theta,Vcmax=Vcmax,PPFD=PAR,Tleaf=Tleaf,Jmax=Jmax,
#                  Rd0=0.68,Q10=2.1,TrefR=22,Rdayfrac=0.7,g1=g1,gsmodel="BBOpti",VPD=VPD,Ca=CO2S,Tcorrect=F)
#   # ret <- ps$Aj - ps$Rd
#   ret <- ps$ALEAF
#   return(ret)
# }

pslow <- function(alpha,PAR,theta,Tleaf,Vcmax,Jmax,Rdayfrac=NULL,g1,gsmodel=NULL,
                  VPD,CO2S,Rd0,delsV,delsJ,EaV,EaJ,EdVC=2e+05,Q10=2) {
  
  ps <- Photosyn(alpha=alpha,theta=theta,Vcmax=Vcmax,PPFD=PAR,Tleaf=Tleaf,Jmax=Jmax,
                 g1=g1,gsmodel="BBOpti",VPD=VPD,Ca=CO2S,Rd0=Rd0,delsC=delsV*10^3,delsJ=delsJ*10^3,EaV=EaV*10^3,EaJ=EaJ*10^3,EdVC=2e+05,Q10=2)
  # ret <- ps$Aj - ps$Rd
  ret <- ps$ALEAF
  return(ret)
}
# 

#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------

# read in-situ photosynthesis data measured on 3rd Feb 2016 (AQ data)

avt.l<-getAQ()

# dat_for_fit<-avt.l


#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
# read Rday (corrected for rday frac =0.7)
rday<-read.csv("Parameters/great_Resp_leaf_area_fh.csv")

avt.l<-merge(avt.l,rday,by="Room")

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------


#--- read g1 estimated from same dataset 

# short_g1<-read.csv("Parameters/great_g1_long_term.csv") # short-term tempearture response

long_g1<-read.csv("Parameters/great_g1_long_term.csv") # long-term tempearture response


#--- read temperature response parameters of Vcmax and Jmax (Ea, delS....)

aci_param<-read.csv("yplantupscale/data/great_aci_t_response.csv")[c(1,3,4,6,7,10,11,13,14)] 


#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------

#--- make a dataframe of physiological parameters

dat_ph<-merge(long_g1,aci_param,by="Room")

#---

#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------


# #-- read original estimates of Vcmax and Jmax (at 25C) estimated from ACi curves
# tres_param<-read.csv("yplantupscale/data/great_aci_t_response.csv")
# tres_param<-subset(tres_param,tres_param$Room %in% c(1,4,6))
# tres_param$JVr<-with(tres_param,Jmax25/Vcmax25)
# 
# # with(tres_param,plot(Tgrowth,JVr))
# # fit a linear model to estimate JVr for other temperature treatments


#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------

#-- need to estimate JV ratio (at each growth temperatures) to get Jmax
#-- use estimates from ACi data to get a relationship

tres_jvr<-read.csv("Parameters/great_jvr_25C.csv") #JV ratio for three growth temperatures 18, 28.5 and 35.5
tres_jvr[2]<-NULL
# lm1<-lm(JVr~Tgrowth,data=tres_jvr)  # fit a linear model 
# abline(summary(lm1)$coefficients[[1]],summary(lm1)$coefficients[[2]])

Tgrowth<-c(22.5,25.0,32.5)
Room<-c(2,3,5)
JVr_rest<-data.frame(Tgrowth,Room)

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
interp_phy <- function(d, Tgrowth,var){
  
  sp <- spline(x=d$Tgrowth, y=d[[var]], xout=c(22.5,25,32.5))
  re <-data.frame(Tgrowth=Tgrowth, var=sp$y)
  return(re)
  
}
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

JVr_rest$JVr<-interp_phy(tres_jvr,Tgrowth=c(22.5,25,32.5),var="JVr")[[2]]

# JVr_rest$JVr<-with(JVr_rest,summary(lm1)$coefficients[[1]]+summary(lm1)$coefficients[[2]]*Tgrowth)

JVr_dat<-rbind(JVr_rest,tres_jvr)
# JVr_dat$JVr_e<-with(JVr_dat,2.32-0.03*Tgrowth)

JVr_dat$Tgrowth<-NULL

with(JVr_dat,plot(Room,JVr))
#-- add other physiological data

dat_ph<-merge(dat_ph,JVr_dat ,by="Room")
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------

#--- add physiological data to the photosynthesis dataset

dat_all<-merge(avt.l,dat_ph,by="Room")

# dat_c1<-subset(dat_all,dat_all$campaign==1 & dat_all$W_treatment=="w") # get data for campaign 1 (2016/02/03)

dat_c1<-dat_all
datlist<-split(dat_c1,dat_c1$Room)

#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------


# estimate Vcmax25 and alpha using nls function
# estimate Vcmax25 and alpha using nls function

alpha<-alpha.se <- theta<-theta.se<-Vcmax25<-Vcmax25.se<-Room <- tair <- c()

for (i in 1:length(unique(dat_c1$Room))) {
  fitalpha <- nls(Photo ~ pslow(alpha, theta,PAR=PARi,Tleaf=Tleaf,Vcmax,Jmax=round(JVr,3)*Vcmax,g1=g1,VPD=VpdL,CO2S=CO2S,Rd0=R_leaf_a*.7,Q10=2,
                                delsV=delsV,delsJ=delsJ,EaV=EaV,EaJ=EaJ),
                  start=list(alpha=0.2,Vcmax=70,theta=.3), data=datlist[[i]], trace=TRUE)
  
  alpha[i]<-summary(fitalpha)$parameters[[1]]
  alpha.se[i]<-summary(fitalpha)$parameters[[4]]
  theta[i]<-summary(fitalpha)$parameters[[3]]
  theta.se[i]<-summary(fitalpha)$parameters[[6]]
  
  Vcmax25[i]<-summary(fitalpha)$parameters[[2]]
  Vcmax25.se[i]<-summary(fitalpha)$parameters[[5]]
  Room[i] <- i
  # tair[i] <- mean(datlist[[i]]$Tair)
}

vj_c1<-data.frame(alpha,alpha.se,theta,theta.se,Vcmax25,Vcmax25.se,Room)
vj_c1$Campaign<-1
vj_c1<-merge(vj_c1,JVr_dat,by="Room")
vj_c1$Jmax25<-vj_c1$Vcmax25*vj_c1$JVr

vj_c1$JVr<-NULL
vj_c1$Date<-as.Date("2016-02-03")
# vj_c1$Jmax25<-ifelse(vj_c1$Room %in% c(5,6),vj_c1$Jmax25*1.3,vj_c1$Jmax25)
# vj_c1$Vcmax25<-ifelse(vj_c1$Room %in% c(5,6),vj_c1$Vcmax25*1.3,vj_c1$Vcmax25)

mean(subset(vj_c1,vj_c1$theta>0)$theta)


# nls(Photo ~ pslow(alpha, theta,PAR=PARi, Tleaf=Tleaf,Vcmax,Jmax= JVr*Vcmax,g1=g1,VPD=VpdL,CO2S=CO2S),
#     start=list(Vcmax=120,alpha=.3,theta=.3), data=datlist[[4]], trace=TRUE)
#---------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------
# mean(subset(vj_c1,vj_c1$theta>0)$theta)
# add g1 to the table

phys_params_all<-merge(dat_ph,vj_c1,by="Room")
phys_params_all$Jmax25.se<-with(phys_params_all,Vcmax25.se) # assume Vcmax SE = Jmax SE


write.csv(phys_params_all,file = "Parameters/great_alpha_vcmax_jmax_25.csv",row.names=FALSE,sep=",") 
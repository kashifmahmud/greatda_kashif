# This sript runs the model equations for parameter shifting from potted seedling to free seedling
# average the daily data
data.attrib.set.daily = summaryBy(Tgrowth+SLA+LA+Leafmass+Stemmass+Rootmass+g1+EaV+delsV+EaJ+delsJ+JVr+alpha+theta+Vcmax25+Jmax25+k+Y+af+as+ar+
                                    Intercept+Slope+self_s.mean ~ Room+Date, data=data.attrib.set, FUN=c(mean), na.rm=TRUE)
names(data.attrib.set.daily) = c("Room","Date","Tgrowth","SLA","LA","LM","SM","RM","g1","EaV","delsV","EaJ","delsJ","JVr","alpha","theta",
                                 "Vcmax25","Jmax25","k","Y","af","as","ar","Intercept","Slope","self_s")


# # adjust the respiration data
# if (data.attrib.set$Room[1] == 1) {
#   c3 = 0.685
# } else {
#   c3 = 0.773705
# }  
  
# plot(data.attrib.set.daily$Date,data.attrib.set.daily$SLA,col='red',type='l')
# lines(data.attrib.set.daily$Date,data.attrib.set.daily$LA/data.attrib.set.daily$LM,col='green')

# Calculating model outputs
Mleaf = Mstem = Mroot = LA = Rm = c()
Mleaf[1] <- data.attrib.set.daily$LM[1]
Mstem[1] <- data.attrib.set.daily$SM[1]
Mroot[1] <- data.attrib.set.daily$RM[1]
LA[1] <- data.attrib.set.daily$LA[1]

k=data.attrib.set.daily$k; Y=data.attrib.set.daily$Y; af=data.attrib.set.daily$af; as=data.attrib.set.daily$as
ar=data.attrib.set.daily$ar; sf=data.attrib.set.daily$sf
Cday = A = GPP = Cstorage = Sleaf = Sstem = Sroot = M = c()


data.attrib.set.date = subset(data.attrib.set, Date %in% as.Date(data.attrib.set.daily$Date[1]))
met.rd.attrib.set.date = subset(met.rd.attrib.set, Date %in% as.Date(data.attrib.set.daily$Date[1]))

A_pred <- Photosyn (VPD=data.attrib.set.date$VPD, Ca=400, Tleaf=data.attrib.set.date$Tair, PPFD=data.attrib.set.date$PAR, gsmodel = "BBOpti", 
                    g1=data.attrib.set.date$g1, alpha = data.attrib.set.date$alpha, theta = data.attrib.set.date$theta, 
                    Jmax = data.attrib.set.date$Jmax25, Vcmax = data.attrib.set.date$Vcmax25, delsC=data.attrib.set.date$delsV*10^3, 
                    delsJ=data.attrib.set.date$delsJ*10^3, EaV=data.attrib.set.date$EaV*10^3, EaJ=data.attrib.set.date$EaJ*10^3, Rd=0)
# need a new dfr with Aleaf and Anet across the day
Aleaf <- A_pred[,c(1:4, 7:11)]
Aleaf_15min <- cbind(Aleaf, data.attrib.set.date[,c("Room","DateTime_hr","Date","Tgrowth")])

Aleaf_15min$Date <- as.Date(Aleaf_15min$Date)
Aleaf_15min$Room <- as.factor(Aleaf_15min$Room)
A[1] <- summaryBy(ALEAF ~ Date+Room, data=Aleaf_15min, FUN=sum, keep.names=TRUE )[1,3] # photosynthetic rate in m-2 s-1 unit

Aleaf_15min$carbon_day <- with(Aleaf_15min, ALEAF*15*60*10^-6*12.0107) # unit conversion from micromole m-2 s-1 to gC m-2 per 15 mins
Cday[1] <- summaryBy(carbon_day ~ Date+Room, data=Aleaf_15min, FUN=sum, keep.names=TRUE )[1,3] # total C gain in gC m-2 d-1 unit

# Multiply with self shading factor
M[1] <- data.attrib.set.date$Slope[1] * LA[1] + data.attrib.set.date$Intercept[1]
GPP[1] <- LA[1] * Cday[1] * M[1] # calculate total daily C gain with self shading in gC d-1 plant-1


##########
# M <- sigma.data$b * LA.data$LA + sigma.data$intercept
# GPP <- LA.data$LA * Cday.data$Cday * M # calculate total daily C gain with self shading
##########

# From Duan's experiment for TNC partitioning to tree organs
# Leaf TNC/Leaf DW =  0.1401421; Stem TNC/Stem DW =  0.0453869; Root TNC/Root DW =  0.02154037
# Sleaf[1] = Mleaf[1] / 0.65 * 0.1401421
# Sstem[1] = Mstem[1] / 0.65 * 0.0453869
# Sroot[1] = Mroot[1] / 0.65 * 0.02154037
Sleaf[1] = Mleaf[1] * 0.1167851
Sstem[1] = Mstem[1] * 0.03782242
Sroot[1] = Mroot[1] * 0.01795031
Cstorage[1] <- Sleaf[1] + Sstem[1] + Sroot[1] 

Cleaf <- Croot <- Cstem <- c()
Cleaf[1] <- Mleaf[1] - Sleaf[1]
Cstem[1] <- Mstem[1] - Sstem[1]
Croot[1] <- Mroot[1] - Sroot[1]

# calculate respiration rates
# met.rd.attrib.set.date[,c("Rdtair_leaf","Rdtair_stem","Rdtair_root")] = 
#   with(met.rd.attrib.set.date, met.rd.attrib.set.date[,c("R25_leaf","R25_stem","R25_root")] * 2.1^((Tair-25)/10)) # unit (gC per gDM per sec)
# met.rd.attrib.set.date[,c("Rdtair_leaf","Rdtair_stem","Rdtair_root")] = 
#   met.rd.attrib.set.date[,c("Rdtair_leaf","Rdtair_stem","Rdtair_root")] * (1/c3) * (15*60) # unit (gC per gC per 15 mins)
# met.rd.attrib.set.date[met.rd.attrib.set.date$period == "Day","Rdtair_leaf"] = 0.7*met.rd.attrib.set.date[met.rd.attrib.set.date$period == "Day","Rdtair_leaf"] # 30% reduction in leaf respiration during day time

met.rd.attrib.set.date$Rdtair_leaf = ifelse(met.rd.attrib.set.date$period=="Day",with(check.rd,Rday_leaf*0.7^((Tair-25)/10)*(LA/Leafmass)), 
                                            with(met.rd.attrib.set.date,R25_leaf*2.1^((Tair-25)/10)))
met.rd.attrib.set.date[,c("Rdtair_stem","Rdtair_root")] = 
  with(met.rd.attrib.set.date, met.rd.attrib.set.date[,c("R25_stem","R25_root")] * 2.1^((Tair-25)/10)) # unit (gC per gC per sec)
met.rd.attrib.set.date[,c("Rdtair_leaf","Rdtair_stem","Rdtair_root")] = 
  met.rd.attrib.set.date[,c("Rdtair_leaf","Rdtair_stem","Rdtair_root")] * (15*60) # unit (gC per gC per 15 mins)

# # Calculate daily mean respiration rates for all tree components by summing all 15-mins data for each day
# Rd.df <- summaryBy(Rdtair_leaf+Rdtair_stem+Rdtair_root ~ Date+Room, data=met.rd.attrib.set.date, FUN=sum, na.rm=TRUE) # Sum of all same day Rd
# 
# Rm[1] = Rd.df$Rdtair_leaf.sum*Mleaf[1] + Rd.df$Rdtair_root.sum*Mroot[1] + Rd.df$Rdtair_stem.sum*Mstem[1]
Rm[1] = sum(met.rd.attrib.set.date$Rdtair_leaf)*Mleaf[1] + sum(met.rd.attrib.set.date$Rdtair_root)*Mroot[1] +
  sum(met.rd.attrib.set.date$Rdtair_stem)*Mstem[1]

# Rm[1] = sum(data.attrib.set.date$R_leaf)*Mleaf[1] + sum(data.attrib.set.date$R_root)*Mroot[1] + 
#   sum(data.attrib.set.date$R_stem)*Mstem[1]
# LA[1] <- LA.data$LA[1]

for (i in 2:nrow(data.attrib.set.daily)) {
  # consider the data file for particular day (i-th day)
  data.attrib.set.date = subset(data.attrib.set, Date %in% as.Date(data.attrib.set.daily$Date[i]))
  # data.attrib.set.date = subset(data.attrib.set, Date %in% as.Date(data.attrib.set$Date[i]))
  met.rd.attrib.set.date = subset(met.rd.attrib.set, Date %in% as.Date(data.attrib.set.daily$Date[i]))
  
  # R_leaf, R_stem and R_root are in gC gC-1 15min-1, so need to convert 15min-1 to d-1
  Cstorage[i] <- Cstorage[i-1] + GPP[i-1] - Rm[i-1] - data.attrib.set.date$k[1]*Cstorage[i-1]
  
  # Cstorage[i] <- Cstorage[i-1] + GPP.data$GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1]) - k[i-1]*Cstorage[i-1]
  Sleaf[i] <- Cstorage[i] * 0.75 # 75% of storage goes to leaf (Duan's experiment)
  Sstem[i] <- Cstorage[i] * 0.16 # 16% of storage goes to stem (Duan's experiment)
  Sroot[i] <- Cstorage[i] * 0.09 # 9% of storage goes to root (Duan's experiment)
  
  Cleaf[i] <- Cleaf[i-1] + data.attrib.set.date$k[1] * Cstorage[i] * data.attrib.set.date$af[1] * (1-data.attrib.set.date$Y[1])
  Cstem[i] <- Cstem[i-1] + data.attrib.set.date$k[1] * Cstorage[i] * data.attrib.set.date$as[1] * (1-data.attrib.set.date$Y[1])
  Croot[i] <- Croot[i-1] + data.attrib.set.date$k[1] * Cstorage[i] * (1-data.attrib.set.date$af[1]-data.attrib.set.date$as[1]) * (1-data.attrib.set.date$Y[1])
  
  Mleaf[i] <- Cleaf[i] + Sleaf[i]
  Mstem[i] <- Cstem[i] + Sstem[i]
  Mroot[i] <- Croot[i] + Sroot[i]
  
  # LA[i] <- data.attrib.set.date$SLA[1] * Mleaf[i] # using SLA from harvest
  # LA[i] <- data.attrib.set.date$LA[1] / data.attrib.set.date$Leafmass[1] * Cleaf[i] # using SLA from modelled LA and LM
  LA[i] <- data.attrib.set.date$LA[1] / data.attrib.set.date$Leafmass[1] * Mleaf[i] # using SLA from modelled LA and LM
  
  # model, should return the aleaf for every 15 minutes. (will retrun Aleaf at 15min interval)
  A_pred <- Photosyn (VPD=data.attrib.set.date$VPD, Ca=400, Tleaf=data.attrib.set.date$Tair, PPFD=data.attrib.set.date$PAR, gsmodel = "BBOpti", 
                      g1=data.attrib.set.date$g1, alpha = data.attrib.set.date$alpha, theta = data.attrib.set.date$theta, 
                      Jmax = data.attrib.set.date$Jmax25, Vcmax = data.attrib.set.date$Vcmax25, delsC=data.attrib.set.date$delsV*10^3, 
                      delsJ=data.attrib.set.date$delsJ*10^3, EaV=data.attrib.set.date$EaV*10^3, EaJ=data.attrib.set.date$EaJ*10^3, Rd=0)
  
  # need a new dfr with Aleaf and Anet across the day
  Aleaf <- A_pred[,c(1:4, 7:11)]
  Aleaf_15min <- cbind(Aleaf, data.attrib.set.date[,c("Room","DateTime_hr","Date","Tgrowth")])
  
  Aleaf_15min$Date <- as.Date(Aleaf_15min$Date)
  Aleaf_15min$Room <- as.factor(Aleaf_15min$Room)
  A[i] <- summaryBy(ALEAF ~ Date+Room, data=Aleaf_15min, FUN=sum, keep.names=TRUE )[1,3] # photosynthetic rate in m-2 s-1 unit
  
  Aleaf_15min$carbon_day <- with(Aleaf_15min, ALEAF*15*60*10^-6*12.0107) # unit conversion from micromole m-2 s-1 to gC m-2 per 15 mins
  Cday[i] <- summaryBy(carbon_day ~ Date+Room, data=Aleaf_15min, FUN=sum, keep.names=TRUE )[1,3] # total GPP in gC m-2 d-1 unit
  
  # Multiply with self shading factor
  M[i] <- data.attrib.set.date$Slope[i] * LA[i] + data.attrib.set.date$Intercept[i]
  GPP[i] <- LA[i] * Cday[i] * M[i] # calculate total daily C gain with self shading in gC d-1 plant-1
  
  # calculate respiration rates
  met.rd.attrib.set.date$Rdtair_leaf = ifelse(met.rd.attrib.set.date$period=="Day",with(check.rd,Rday_leaf*0.7^((Tair-25)/10)*(LA/Leafmass)), 
                                              with(met.rd.attrib.set.date,R25_leaf*2.1^((Tair-25)/10)))
  met.rd.attrib.set.date[,c("Rdtair_stem","Rdtair_root")] = 
    with(met.rd.attrib.set.date, met.rd.attrib.set.date[,c("R25_stem","R25_root")] * 2.1^((Tair-25)/10)) # unit (gC per gC per sec)
  met.rd.attrib.set.date[,c("Rdtair_leaf","Rdtair_stem","Rdtair_root")] = 
    met.rd.attrib.set.date[,c("Rdtair_leaf","Rdtair_stem","Rdtair_root")] * (15*60) # unit (gC per gC per 15 mins)
  
  # met.rd.attrib.set.date[,c("Rdtair_leaf","Rdtair_stem","Rdtair_root")] = 
  #   with(met.rd.attrib.set.date, met.rd.attrib.set.date[,c("R25_leaf","R25_stem","R25_root")] * 2.1^((Tair-25)/10)) # unit (gC per gDM per sec)
  # met.rd.attrib.set.date[,c("Rdtair_leaf","Rdtair_stem","Rdtair_root")] = 
  # met.rd.attrib.set.date[,c("Rdtair_leaf","Rdtair_stem","Rdtair_root")] * (1/c3) * (15*60) # unit (gC per gC per 15 mins)
  # met.rd.attrib.set.date[met.rd.attrib.set.date$period == "Day","Rdtair_leaf"] = 0.7*met.rd.attrib.set.date[met.rd.attrib.set.date$period == "Day","Rdtair_leaf"] # 30% reduction in leaf respiration during day time
  
  # # Calculate daily mean respiration rates for all tree components by summing all 15-mins data for each day
  # Rd.df = c()
  # Rd.df <- summaryBy(Rdtair_leaf+Rdtair_stem+Rdtair_root ~ Date+Room, data=met.rd.attrib.set.date, FUN=sum, na.rm=TRUE) # Sum of all same day Rd
  # 
  # Rm[i] = Rd.df$Rdtair_leaf.sum*Mleaf[i] + Rd.df$Rdtair_root.sum*Mroot[i] + Rd.df$Rdtair_stem.sum*Mstem[i]
  Rm[i] = sum(met.rd.attrib.set.date$Rdtair_leaf)*Mleaf[i] + sum(met.rd.attrib.set.date$Rdtair_root)*Mroot[i] +
    sum(met.rd.attrib.set.date$Rdtair_stem)*Mstem[i]
  # Rm[i] = sum(data.attrib.set.date$R_leaf)*Mleaf[i] + sum(data.attrib.set.date$R_root)*Mroot[i] + 
  #     sum(data.attrib.set.date$R_stem)*Mstem[i]
}

# if (q == 0) {
#   c2 = (output.attrib.set$Mleaf.modelled[i]+output.attrib.set$Mwood.modelled[i]+output.attrib.set$Mroot.modelled[i]) / 
#     (Mleaf[i]+Mstem[i]+Mroot[i])
# }
# if (q > 0) {
#   if (data.attrib.set$Room[1] == 1) {
#     c2 = 1.001311
#   } else {
#     c2 = 1.259885
#   }
# }
# adjust the biomass data
if (data.attrib.set$Room[1] == 1 && q > 0) {
  c2 = c2 + ((6.654/5.29)-1)/12
} else if (data.attrib.set$Room[1] == 4 && q == 0) {
  c2 = 6.654/5.68
} else if (data.attrib.set$Room[1] == 4 && q == 1) {
  c2 = c2
} else if (data.attrib.set$Room[1] == 4 && q == 2) {
  c2 = c2
} else if (data.attrib.set$Room[1] == 4 && q > 2) {
  c2 = c2 + ((3.2446/2.32)-(6.654/5.68))/10
}  else {
  c2 = 1
}
# if (data.attrib.set$Room[1] == 1 && q > 0) {
#   c2 = c2 + ((6.654/5.29)-1)/12
# } else if (data.attrib.set$Room[1] == 4 && q == 0) {
#   c2 = 6.654/5.68
# } else if (data.attrib.set$Room[1] == 4 && q > 0) {
#   c2 = 3.2446/2.32
# }  else {
#   c2 = 1
# }
output.final = data.frame(c2*Cstorage,c2*Mleaf,c2*Mstem,c2*Mroot,c2*Sleaf,c2*A)
# output.final = data.frame(Cstorage,Mleaf,Mstem,Mroot,Sleaf,A)


# check on first day GPP for room 1 from Dushan's analysis
gpp_dushan_day1 = sum(subset(data.attrib.set, Date %in% as.Date(data.attrib.set.daily$Date[1]))$GPP)*(data.attrib.set.daily$Slope[1]*data.attrib.set.daily$LA[1]+data.attrib.set.daily$Intercept[1])
gpp_dushan_day48 = sum(subset(data.attrib.set, Date %in% as.Date(data.attrib.set.daily$Date[48]))$GPP)
gpp_tot_dushan = sum(data.attrib.set$GPP)
gpp_tot_attrib = sum(GPP)

# Plant Carbon pools for various parameter sensitivity
output.final$Date = data.attrib.set.daily$Date
names(output.final) = c("Cstorage","Mleaf","Mstem","Mroot","Sleaf","A","Date")
melted.output = reshape2::melt(output.final[,c("Mleaf","Mstem","Mroot","Cstorage","Sleaf","A","Date")], id.vars="Date")
melted.Cstorage = output.final[,c("Cstorage","Date")]
melted.output$Date = as.Date(melted.output$Date)
melted.Cstorage$Date = as.Date(melted.Cstorage$Date)
melted.output$Case = as.factor(q)
melted.Cstorage$Case = as.factor(q)

# Storing the summary of data, outputs, Cstorage, parameters
if (q == 0) {
  shift.output = melted.output
  shift.Cstorage = melted.Cstorage
}
if (q > 0) {
  shift.output = rbind(shift.output,melted.output)
  shift.Cstorage = rbind(shift.Cstorage,melted.Cstorage)
}

# This code performs sensitivity analysis with parameter shifting from Room 4 (25.8 degree C) to Room 6 (35.5 degree C)
# and make Figure 
#-------------------------------------------------------------------------------------
# Load the custom analysis and plotting functions that do all of the actual work	
source("R/functions_great.R")	
source("R/functions_great_CBM.R")	

# Read tree attributes data including Met data, Temperature dependant variables, Modelled Parameters from Dushan's analysis 
data.attrib <- read.csv("parameters/data_for_attribution_analysis_v2.csv")
data.attrib$Tgrowth[data.attrib$Tgrowth == 25.8] = 28.5 #correction to temperature from 25.8 to 25.5 

data.attrib$Date = as.Date(data.attrib$Date, format="%d/%m/%Y")
# data.attrib$Date = as.Date(data.attrib$Date)
data.attrib$biomass = data.attrib$Leafmass + data.attrib$Stemmass + data.attrib$Rootmass
data.attrib$Intercept = mean(data.attrib$Intercept)
data.attrib$Slope = mean(data.attrib$Slope)
# data.attrib$SLA = mean(data.attrib$SLA)

# check the self-shading factors
png(file = "output/1.self-shading_comparison.png", units="px", width=500, height=500, res=100)
plot(unique(data.attrib$Room), unique(data.attrib$self_s.mean),xlab="Rooms",ylab="Self-shading factor")
dev.off()

##################------------------------------
# Look for the variation in SLA in Room 1
png(file = "output/1.SLA_comparison.png", units="px", width=1500, height=2500, res=300)
par(mfrow = c(3,1))
data.attrib.set = subset(data.attrib, Room %in% c(1)) # Consider the free seedling to test the parameter sensitivity

data.attrib.set.daily = summaryBy(Tgrowth+SLA+LA+Leafmass+Stemmass+Rootmass+g1+EaV+delsV+EaJ+delsJ+JVr+alpha+theta+Vcmax25+Jmax25+k+Y+af+as+ar+
                                    Intercept+Slope+self_s.mean ~ Room+Date, data=data.attrib.set, FUN=c(mean), na.rm=TRUE)
names(data.attrib.set.daily) = c("Room","Date","Tgrowth","SLA","LA","LM","SM","RM","g1","EaV","delsV","EaJ","delsJ","JVr","alpha","theta",
                                 "Vcmax25","Jmax25","k","Y","af","as","ar","Intercept","Slope","self_s")

plot(data.attrib.set.daily$Date,data.attrib.set.daily$SLA,col='red',type='l',main='Room 1',xlab='',ylab='SLA (m2 gC-1)')
lines(data.attrib.set.daily$Date,data.attrib.set.daily$LA/data.attrib.set.daily$LM,col='green')
legend('topleft', c("SLA from harvest", "Modelled SLA from allometry (LA/LM)"), lty=1, col=c('red','green'), bty='n', cex=0.75)

# Look for the variation in SLA in Room 4
data.attrib.set = subset(data.attrib, Room %in% c(4)) # Consider the free seedling to test the parameter sensitivity

data.attrib.set.daily = summaryBy(Tgrowth+SLA+LA+Leafmass+Stemmass+Rootmass+g1+EaV+delsV+EaJ+delsJ+JVr+alpha+theta+Vcmax25+Jmax25+k+Y+af+as+ar+
                                    Intercept+Slope+self_s.mean ~ Room+Date, data=data.attrib.set, FUN=c(mean), na.rm=TRUE)
names(data.attrib.set.daily) = c("Room","Date","Tgrowth","SLA","LA","LM","SM","RM","g1","EaV","delsV","EaJ","delsJ","JVr","alpha","theta",
                                 "Vcmax25","Jmax25","k","Y","af","as","ar","Intercept","Slope","self_s")

plot(data.attrib.set.daily$Date,data.attrib.set.daily$SLA,col='red',type='l',main='Room 4',xlab='',ylab='SLA (m2 gC-1)')
lines(data.attrib.set.daily$Date,data.attrib.set.daily$LA/data.attrib.set.daily$LM,col='green')
legend('topleft', c("SLA from harvest", "Modelled SLA from allometry (LA/LM)"), lty=1, col=c('red','green'), bty='n', cex=0.75)

# Look for the variation in SLA in Room 6
data.attrib.set = subset(data.attrib, Room %in% c(6)) # Consider the free seedling to test the parameter sensitivity

data.attrib.set.daily = summaryBy(Tgrowth+SLA+LA+Leafmass+Stemmass+Rootmass+g1+EaV+delsV+EaJ+delsJ+JVr+alpha+theta+Vcmax25+Jmax25+k+Y+af+as+ar+
                                    Intercept+Slope+self_s.mean ~ Room+Date, data=data.attrib.set, FUN=c(mean), na.rm=TRUE)
names(data.attrib.set.daily) = c("Room","Date","Tgrowth","SLA","LA","LM","SM","RM","g1","EaV","delsV","EaJ","delsJ","JVr","alpha","theta",
                                 "Vcmax25","Jmax25","k","Y","af","as","ar","Intercept","Slope","self_s")

plot(data.attrib.set.daily$Date,data.attrib.set.daily$SLA,col='red',type='l',main='Room 6',xlab='',ylab='SLA (m2 gC-1)')
lines(data.attrib.set.daily$Date,data.attrib.set.daily$LA/data.attrib.set.daily$LM,col='green')
legend('topleft', c("SLA from harvest", "Modelled SLA from allometry (LA/LM)"), lty=1, col=c('red','green'), bty='n', cex=0.75)

dev.off()

##################------------------------------
# take only Room 4 and 6 for attribution analysis
data.attrib = subset(data.attrib, Room %in% c(4,6))


##################------------------------------
# Consider everything (Cday, LA, Rd, Biomass, parameters) for Room 4
q=0 # Case 0
data.attrib.set = subset(data.attrib, Room %in% c(4)) # Start with Room 4 to test the parameter sensitivity
# data.attrib.set$SLA = (subset(data.attrib, Room %in% c(6)))$SLA # Consider SLA from Room 6 to match the attribution
# data.attrib.set$LA = (subset(data.attrib, Room %in% c(6)))$LA # Consider LA from Room 6 to match the attribution
# data.attrib.set$Leafmass = (subset(data.attrib, Room %in% c(6)))$Leafmass # Consider LM from Room 6 to match the attribution

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# Meteorological data attribution
q=1 # Case 1
# Modify the Vcmax and Jmax for Room 6
data.attrib.set$VPD = (subset(data.attrib, Room %in% c(6)))$VPD
# data.attrib.set$Tair = (subset(data.attrib, Room %in% c(6)))$Tair
# data.attrib.set$PAR = (subset(data.attrib, Room %in% c(6)))$PAR

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# Meteorological data attribution
q=2 # Case 2
# Modify the Vcmax and Jmax for Room 6
# data.attrib.set$VPD = (subset(data.attrib, Room %in% c(6)))$VPD
data.attrib.set$Tair = (subset(data.attrib, Room %in% c(6)))$Tair
# data.attrib.set$PAR = (subset(data.attrib, Room %in% c(6)))$PAR

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# Vcmax and Jmax attribution
q=3 # Case 3
# Modify the Vcmax and Jmax for Room 6
data.attrib.set$Vcmax25 = (subset(data.attrib, Room %in% c(6)))$Vcmax25
data.attrib.set$Jmax25 = (subset(data.attrib, Room %in% c(6)))$Jmax25

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# Temperature response parameters of Vcmax and Jmax (Activation energy (Ea), entropy (∆S))
q=4 # Case 4
# Modify the Ea and ∆S for Room 6
data.attrib.set$EaV = (subset(data.attrib, Room %in% c(6)))$EaV
data.attrib.set$delsC = (subset(data.attrib, Room %in% c(6)))$delsC
data.attrib.set$EaJ = (subset(data.attrib, Room %in% c(6)))$EaJ
data.attrib.set$delsJ = (subset(data.attrib, Room %in% c(6)))$delsJ

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# g1
q=5 # Case 5
# Modify the g1 for Room 6
data.attrib.set$g1 = (subset(data.attrib, Room %in% c(6)))$g1

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# alpha and theta
q=6 # Case 6
# Modify the alpha for Room 6
data.attrib.set$alpha = (subset(data.attrib, Room %in% c(6)))$alpha
data.attrib.set$theta = (subset(data.attrib, Room %in% c(6)))$theta

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# Rm (leaf, stem and root)
q=7 # Case 7
# Modify the R_leaf, R_stem and R_root for Room 6
data.attrib.set$R_leaf = (subset(data.attrib, Room %in% c(6)))$R_leaf
data.attrib.set$R_stem = (subset(data.attrib, Room %in% c(6)))$R_stem
data.attrib.set$R_root = (subset(data.attrib, Room %in% c(6)))$R_root

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# SLA
q=8 # Case 8
# Modify the SLA for Room 6
data.attrib.set$SLA = (subset(data.attrib, Room %in% c(6)))$SLA # Consider SLA from Room 6 to match the attribution
data.attrib.set$LA = (subset(data.attrib, Room %in% c(6)))$LA # Consider LA from Room 6 to match the attribution
data.attrib.set$Leafmass = (subset(data.attrib, Room %in% c(6)))$Leafmass # Consider LM from Room 6 to match the attribution

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# af, as, ar
q=9 # Case 9
# Modify the af, as and ar for Room 6
data.attrib.set$af = (subset(data.attrib, Room %in% c(6)))$af
data.attrib.set$as = (subset(data.attrib, Room %in% c(6)))$as
data.attrib.set$ar = (subset(data.attrib, Room %in% c(6)))$ar

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# Y
q=10 # Case 10
# Modify the af, as and ar for Room 6
data.attrib.set$Y = (subset(data.attrib, Room %in% c(6)))$Y

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# k
q=11 # Case 11
# Modify the af, as and ar for Room 6
data.attrib.set$k = (subset(data.attrib, Room %in% c(6)))$k

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

# ##################------------------------------
# # self shading factors
# q=11 # Case 11
# # Modify the Slope and Intercept for Room 6
# data.attrib.set$Slope = (subset(data.attrib, Room %in% c(6)))$Slope
# data.attrib.set$Intercept = (subset(data.attrib, Room %in% c(6)))$Intercept
# 
# # This sript runs the model equations for parameter shifting from potted seedling to free seedling
# source("R/CBM_model_shift_great.R")


############ Summarize the C pools plots
# Load the custom analysis and plotting functions that do all of the actual work	
source("R/functions_great.R")	

plot.shift = list() 
font.size = 7
title = as.character(c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"))
cbPalette1 = c("skyblue", "orange")

data.attrib$Date = as.Date(data.attrib$Date)
data.attrib.daily = summaryBy(VPD+Tair+PAR+Tgrowth+SLA+LA+Leafmass+Stemmass+Rootmass+g1+EaV+delsV+EaJ+delsJ+JVr+alpha+theta+Vcmax25+Jmax25+k+Y+af+as+ar+
                                Intercept+Slope+self_s.mean ~ Room+Date, data=data.attrib, FUN=c(mean), na.rm=TRUE)
data.attrib.daily.resp = summaryBy(R_leaf+R_stem+R_root ~ Room+Date, data=data.attrib, FUN=c(sum), na.rm=TRUE)
data.attrib.daily = merge(data.attrib.daily, data.attrib.daily.resp, by=c("Room", "Date"))
names(data.attrib.daily) = c("Room","Date","VPD","Tair","PAR","Tgrowth","SLA","LA","LM","SM","RM","g1","EaV","delsV","EaJ","delsJ","JVr","alpha","theta",
                             "Vcmax25","Jmax25","k","Y","af","as","ar","Intercept","Slope","self_s","R_leaf","R_stem","R_root")

######## Plot Meteorological data
plot.shift[[1]] = plot.VPD.v2(data.attrib.daily, 1)
plot.shift[[2]] = plot.Tair(data.attrib.daily, 1)
# plot.shift[[3]] = plot.PAR(data.attrib.daily, 3)

######## Plot both Vcmax25 and Jmax25
plot.shift[[3]] = plot.Vcmax.Jmax.v2(data.attrib.daily, 1)
# plot.shift[[3]] = plot.Jmax(data.attrib.daily, 5)

######## Plot Activation energy (Ea) and entropy (∆S)
# get Vcmax temperature response
tleaf<-seq(15,45,1)
data.attrib.daily.rm4 = subset(data.attrib.daily, Room == 4) 
room4 = fit_parr_vs(Ts=tleaf,Ea=mean(data.attrib.daily.rm4$EaV),delS=mean(data.attrib.daily.rm4$delsV))
Vcmax.response.rm4 = data.frame(tleaf=tleaf, Room = 4, response = room4)
data.attrib.daily.rm6 = subset(data.attrib.daily, Room == 6) 
room6 = fit_parr_vs(Ts=tleaf,Ea=mean(data.attrib.daily.rm6$EaV),delS=mean(data.attrib.daily.rm6$delsV))
Vcmax.response.rm6 = data.frame(tleaf=tleaf, Room = 6, response = room6)
Vcmax.response = rbind(Vcmax.response.rm4,Vcmax.response.rm6)

# plot.shift[[4]] = plot.Ea(data.attrib.daily, 1)
plot.shift[[4]] = plot.Ea(Vcmax.response, 1)

# get Jmax temperature response
data.attrib.daily.rm4 = subset(data.attrib.daily, Room == 4) 
room4 = fit_parr_vs(Ts=tleaf,Ea=mean(data.attrib.daily.rm4$EaJ),delS=mean(data.attrib.daily.rm4$delsJ))
Jmax.response.rm4 = data.frame(tleaf=tleaf, Room = 4, response = room4)
data.attrib.daily.rm6 = subset(data.attrib.daily, Room == 6) 
room6 = fit_parr_vs(Ts=tleaf,Ea=mean(data.attrib.daily.rm6$EaJ),delS=mean(data.attrib.daily.rm6$delsJ))
Jmax.response.rm6 = data.frame(tleaf=tleaf, Room = 6, response = room6)
Jmax.response = rbind(Jmax.response.rm4,Jmax.response.rm6)

# plot.shift[[5]] = plot.dels(data.attrib.daily, 1)
plot.shift[[5]] = plot.dels(Jmax.response, 1)

######## Plot g1
plot.shift[[6]] = plot.g1(data.attrib.daily, 1)

######## Plot alpha and theta
PAR<-seq(-100,2000,10)
data.attrib.daily.rm4 = subset(data.attrib.daily, Room == 4) 
room4 = curve.nlslrc(data.attrib.daily.rm4$alpha, data.attrib.daily.rm4$theta, mean(data.attrib.daily.rm4$Jmax25))
light.response.rm4 = data.frame(PAR=PAR, Room = 4, response = room4$J)
data.attrib.daily.rm6 = subset(data.attrib.daily, Room == 6) 
room6 = curve.nlslrc(data.attrib.daily.rm6$alpha, data.attrib.daily.rm6$theta, mean(data.attrib.daily.rm6$Jmax25))
light.response.rm6 = data.frame(PAR=PAR, Room = 6, response = room6$J)
light.response = rbind(light.response.rm4,light.response.rm6)

plot.shift[[7]] = plot.alpha.theta(light.response, 1)

######## Plot Respiration rates
plot.shift[[8]] = plot.Rd(data.attrib.daily, 1)

######## Plot SLA
sla.harvest.all = read.csv("processed_data/sla.harvest.all.csv")
sla.harvest.all = subset(sla.harvest.all, Room %in% as.factor(c("4","6")))
sla.harvest.all$Date = as.Date(sla.harvest.all$Date)

plot.shift[[9]] = plot.SLA(data.attrib.daily, sla.harvest.all, 1)

# Plot individual modelled parameters ("k","Y","af","sf") against "volume"
plot.shift[[10]] = plot.allocation.fractions(data.attrib.daily, 1)

plot.shift[[11]] = plot.Y(data.attrib.daily, 1)

plot.shift[[12]] = plot.k(data.attrib.daily, 1)

png("output/Figure_parameter_shifting_rm4-6_params.png", units="px", width=3000, height=1500, res=250)
lay <- rbind(c(1,2,3),c(4,5,6),c(7,8,9),c(10,11,12))
grid.arrange(grobs = plot.shift, layout_matrix = lay)
dev.off()

# # This sript plots the C pools for various test cases with parameter shifted from potted seedling to free seedling
# source("R/generate_figures_param_shift.R")
font.size = 12
# cbPalette = c("gray", "orange", "skyblue", "green", "black", "yellow3", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "mediumorchid", "lightpink")
cbPalette = c("gray", "orange", "skyblue", "green", "black", "yellow3", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "mediumorchid", "olivedrab1")

shift.output.Mleaf = subset(shift.output,(variable %in% "Mleaf"))
plot.shift[[13]] = plot.Mleaf.rm4to6(shift.output.Mleaf)

shift.output.Mstem = subset(shift.output,(variable %in% "Mstem"))
plot.shift[[14]] = plot.Mstem.rm4to6(shift.output.Mstem)

shift.output.Mroot = subset(shift.output,(variable %in% "Mroot"))
plot.shift[[15]] = plot.Mroot.rm4to6(shift.output.Mroot)

png("output/Figure_parameter_shifting_rm4-6_linegraphs.png", units="px", width=3000, height=3500, res=300)
lay <- rbind(c(1,2,13,13),c(3,4,13,13),c(5,6,14,14),c(7,8,14,14),c(9,10,15,15),c(11,12,15,15))
grid.arrange(grobs = plot.shift, layout_matrix = lay)
dev.off()

shift.output.dcast = dcast(shift.output, Date+Case ~ variable)
shift.output.dcast$biomass = shift.output.dcast$Mleaf + shift.output.dcast$Mstem + shift.output.dcast$Mroot
# plot.shift[[16]] = plot.biomass.rm1to4(shift.output.dcast)
plot.shift[[16]] = plot.biomass.logscale.rm1to4(shift.output.dcast)

shift.output.bar = subset(shift.output.dcast, Date %in% as.Date("2016-02-24"))
keeps = c("Date","Case","biomass")
shift.output.bar = shift.output.bar[ , keeps, drop = FALSE]
# shift.output.bar[,c("Date","Case","Cstorage","Sleaf")] = NULL
row.names(shift.output.bar) <- as.numeric(paste(shift.output.bar$Case))
# shift.output.bar[,c("Date","Case")] = NULL
plot.shift[[17]] = plot.biomass.barplot.rm4to6(shift.output.bar)
plot2 = plot.shift[[17]] +
  annotate("text", x = max(as.numeric(shift.output.bar$Case)), y = max(shift.output.bar$biomass), size = font.size-7, label = paste("(b)"))

png("output/Figure_parameter_shifting_rm4-6_barchart.png", units="px", width=1000, height=800, res=200)
plot.shift[[17]]
dev.off()

png("output/Figure_parameter_shifting_rm4-6.png", units="px", width=3000, height=4000, res=300)
lay <- rbind(c(1,2,13,13),c(3,4,13,13),c(5,6,14,14),c(7,8,14,14),c(9,10,15,15),c(11,12,15,15),c(16,16,17,17),c(16,16,17,17))
grid.arrange(grobs = plot.shift, layout_matrix = lay)
dev.off()

png("output/Figure_parameter_shifting_barchart.png", units="px", width=1500, height=800, res=150)
grid.arrange(plot1, plot2, ncol=2)
dev.off()

#-------------------------------------------------------------------------------------


# #-------------------------------------------------------------------------------------
# # Quantify the changes of final biomasses
# mass.quantify = subset(shift.output, Date %in% as.Date("2013-05-21"))
# keeps = c("variable", "value", "Case")
# mass.quantify = mass.quantify[ , keeps, drop = FALSE]
# mass.quantify = dcast( mass.quantify , variable ~ Case )
# mass.quantify = mass.quantify[1:3,]
# 
# mass.quantify = mass.quantify[,-1]
# mass.quantify = data.frame(t(mass.quantify))
# # mass.quantify = mass.quantify[-1,]
# rownames(mass.quantify) = c("5L", "Cday", "Rd", "allocations", "Y", "sf", "k")
# names(mass.quantify) = c("Ct,f", "Ct,w", "Ct,r")
# 
# # for (i in 2:7) {
# #   mass.quantify[i,] = mass.quantify[i,] - mass.quantify[1,]
# # }
# # mass.quantify[8,] = colSums(mass.quantify)
# mass.quantify[,4] = rowSums(mass.quantify)
# # rownames(mass.quantify)[8] = "Total"
# names(mass.quantify)[4] = "Ct"
# mass.quantify[,5] = mass.quantify[,4]
# for (i in 7:2) {
#    mass.quantify[i,5] = mass.quantify[i,5] - mass.quantify[(i-1),5]
# }
# 
# # mass.quantify$change.Cday = (mass.quantify[,3] - mass.quantify[,2])
# # mass.quantify$change.Rd = (mass.quantify[,4] - mass.quantify[,2])
# # mass.quantify$change.alloc = (mass.quantify[,5] - mass.quantify[,2])
# # mass.quantify$change.Y = (mass.quantify[,6] - mass.quantify[,2])
# # mass.quantify$change.sf = (mass.quantify[,7] - mass.quantify[,2])
# # mass.quantify$change.k = (mass.quantify[,8] - mass.quantify[,2])
# # 
# # mass.quantify$change.Cday = (mass.quantify[,3] - mass.quantify[,2]) / mass.quantify[,2]
# # mass.quantify$change.Rd = (mass.quantify[,4] - mass.quantify[,3]) / mass.quantify[,3]
# # mass.quantify$change.alloc = (mass.quantify[,5] - mass.quantify[,4]) / mass.quantify[,4]
# # mass.quantify$change.Y = (mass.quantify[,6] - mass.quantify[,5]) / mass.quantify[,5]
# # mass.quantify$change.sf = (mass.quantify[,7] - mass.quantify[,6]) / mass.quantify[,6]
# # mass.quantify$change.k = (mass.quantify[,8] - mass.quantify[,7]) / mass.quantify[,7]
# # mass.quantify$total = rowSums(mass.quantify[,c(2,9:14)])
# # mass.quantify[4,c(2:8)] = colSums(mass.quantify[,c(2:8)])
# 
# # mass.quantify = mass.quantify[,-c(2:8)]
# write.csv(mass.quantify, file = "output/final_mass_changes.csv", row.names = FALSE)
# #-------------------------------------------------------------------------------------
# 
# 
# #-------------------------------------------------------------------------------------
# #- Make figure for AGU presentation
# source("R/Parameter_shifting_AGU_presentation.R")
# 
# #-------------------------------------------------------------------------------------


# This code performs sensitivity analysis with parameter shifting from Room 1 (18.5 degree C) to Room 4 (25.8 degree C)
# and make Figure 
#-------------------------------------------------------------------------------------
# Load the custom analysis and plotting functions that do all of the actual work	
source("R/functions_great.R")	
source("R/functions_great_CBM.R")	

# Read tree attributes data including Met data, Temperature dependant variables, Modelled Parameters from Dushan's analysis 
# data.attrib <- read.csv("Parameters/data_for_attribution_analysis_v2.csv")
# data.attrib$Tgrowth[data.attrib$Tgrowth == 25.8] = 28.5 #correction to temperature from 25.8 to 28.5 
# data.attrib$Date = as.Date(data.attrib$Date, format="%d/%m/%Y")

data.attrib <- read.csv("Parameters/data.attrib.merged.csv")
data.attrib$Date = as.Date(data.attrib$Date)
data.attrib$DateTime_hr = as.POSIXct(as.character(data.attrib$DateTime_hr), format="%Y-%m-%d %H:%M:%S")

data.attrib$biomass = data.attrib$Leafmass + data.attrib$Stemmass + data.attrib$Rootmass
data.attrib$Intercept = mean(data.attrib$Intercept)
data.attrib$Slope = mean(data.attrib$Slope)
# data.attrib$SLA = mean(data.attrib$SLA)

# Read Rd at 25C from Dushan's analysis (Unit = gC per gC per sec)
Rd.attrib <- read.csv("parameters/Rd_data_for_attribution_analysis.csv")
names(Rd.attrib) = c("Room","R25_leaf","R25_root","R25_stem","R25_leaf_SE","R25_root_SE","R25_stem_SE")

# Read Rday at 25C from Dushan's analysis (Unit = mumol m-2 s-1)
Rday.attrib <- read.csv("parameters/Rday_aci_fits.csv")
keeps = c("Room", "Rd.mean")
Rday.attrib = Rday.attrib[ , keeps, drop = FALSE]
names(Rday.attrib) = c("Room","Rday_leaf")
Rday.attrib$Rday_leaf = Rday.attrib$Rday_leaf*10^-6*12.0107 # Unit converted to gC per gC per sec
Rd.attrib = merge(Rd.attrib,Rday.attrib, by="Room")
Rd.attrib[,c("R25_leaf","Rday_leaf")] = Rd.attrib[,c("R25_leaf","Rday_leaf")]
met.rd.attrib = merge(data.attrib[,c("Date","Room","Tair","DateTime_hr","period","LA","Leafmass")],Rd.attrib, by="Room")

# check respiration rates
check.rd = merge(data.attrib[,c("Date","Room","Tair","DateTime_hr","period","R_leaf","R_stem","R_root","LA","Leafmass")],Rd.attrib, by="Room")
check.rd$Rdtair_leaf = ifelse(check.rd$period=="Day",with(check.rd,Rday_leaf*0.7^((Tair-25)/10)*(LA/Leafmass)), 
                              with(check.rd,R25_leaf*2.1^((Tair-25)/10)))
check.rd[,c("Rdtair_stem","Rdtair_root")] = 
  with(check.rd, check.rd[,c("R25_stem","R25_root")] * 2.1^((Tair-25)/10)) # unit (gC per gC per sec)
# check.rd[,c("Rdtair_leaf","Rdtair_stem","Rdtair_root")] = 
#   check.rd[,c("Rdtair_leaf","Rdtair_stem","Rdtair_root")] * (1/c1) * (15*60) # unit (gC per gC per 15 mins)
check.rd[,c("Rdtair_leaf","Rdtair_stem","Rdtair_root")] = 
  check.rd[,c("Rdtair_leaf","Rdtair_stem","Rdtair_root")] * (15*60) # unit (gC per gC per 15 mins)
# check.rd[check.rd$period == "Day","R_leaf"] = 0.7*check.rd[check.rd$period == "Day","R_leaf"] # 30% reduction in leaf respiration during day time
Rd.df <- summaryBy(R_leaf+R_stem+R_root+Rdtair_leaf+Rdtair_stem+Rdtair_root ~ Date+Room, data=check.rd, FUN=sum, na.rm=TRUE) # Sum of all same day Rd

write.csv(Rd.df, "parameters/Rd.df.csv", row.names=FALSE)


##################------------------------------
# Read parameters
if (length (result[[1]]) == 1) {
  no.param.par.var = result[[1]]
  summary.param = result[[2]]
  summary.data = result[[3]]
  summary.output = result[[4]]
  summary.error = result[[5]]
} else {
  no.param.par.var = result[[1]][[1]]
  summary.param = result[[1]][[2]]
  summary.data = result[[1]][[3]]
  summary.output = result[[1]][[4]]
  summary.error = result[[1]][[5]]
}
keeps = c("Date", "variable", "value", "treatment")
output.set = summary.output[ , keeps, drop = FALSE]
output.attrib = reshape2::dcast( output.set , Date+treatment ~ variable, value.var = "value")
names(output.attrib)[2] = "Room"

##################------------------------------
# check the self-shading factors
png(file = "output/1.self-shading_comparison.png", units="px", width=500, height=500, res=100)
plot(unique(data.attrib$Room), unique(data.attrib$self_s.mean),xlab="Rooms",ylab="Self-shading factor")
dev.off()

# check LA
# cbPalette2 = c("gray", "yellow3", "skyblue", "green3", "#0072B2", "#D55E00", "orange")
cbPalette2 = c("gray","yellow","cyan","blue","magenta","red")
# png(file = "output/1.LA_comparison.png", units="px", width=500, height=500, res=100)
p1 = ggplot() +
  geom_line(data = data.attrib, aes(x = Date, y = LA, group = as.factor(Room), colour=as.factor(Room))) +
  geom_point(size=2) +
  ylab(expression(LA~(m^{2}~plant^{-1}))) + 
  # annotate("text", x = min(shift.output.Mleaf$Date), y = max(shift.output.Mleaf$value), size = font.size, label = paste(title[7])) +
  scale_colour_manual(breaks=c("1","2","3","4","5","6"),labels=c("Room 1","Room 2","Room 3","Room 4","Room 5","Room 6"),values=cbPalette2) +
  theme_bw() +
  theme(legend.position = c(0.25,0.67),legend.text.align = 0) +
  theme(legend.title = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(legend.key.height=unit(1,"line")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

p2 = ggplot() +
  geom_line(data = data.attrib, aes(x = Date, y = Leafmass, group = as.factor(Room), colour=as.factor(Room))) +
  geom_point(size=2) +
  ylab(expression(LM~(g~C~plant^{-1}))) + 
  # annotate("text", x = min(shift.output.Mleaf$Date), y = max(shift.output.Mleaf$value), size = font.size, label = paste(title[7])) +
  scale_colour_manual(breaks=c("1","2","3","4","5","6"),labels=c("Room 1","Room 2","Room 3","Room 4","Room 5","Room 6"),values=cbPalette2) +
  theme_bw() +
  theme(legend.position = c(0.25,0.67),legend.text.align = 0) +
  theme(legend.title = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(legend.key.height=unit(1,"line")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(color=FALSE)

p3 = ggplot() +
  geom_line(data = data.attrib, aes(x = Date, y = LA/Leafmass, group = as.factor(Room), colour=as.factor(Room))) +
  geom_point(size=2) +
  ylab(expression(SLA~(m^{2}~g~C^{-1}))) + 
  # annotate("text", x = min(shift.output.Mleaf$Date), y = max(shift.output.Mleaf$value), size = font.size, label = paste(title[7])) +
  scale_colour_manual(breaks=c("1","2","3","4","5","6"),labels=c("Room 1","Room 2","Room 3","Room 4","Room 5","Room 6"),values=cbPalette2) +
  theme_bw() +
  theme(legend.position = c(0.25,0.67),legend.text.align = 0) +
  theme(legend.title = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=font.size)) +
  theme(legend.key.height=unit(1,"line")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = font.size, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + guides(color=FALSE)

png("output/1.Figure_modelled_sla_great.png", units="px", width=2000, height=700, res=150)
# print (do.call(grid.arrange,  p))
grid.arrange(p1, p2, p3, ncol=3)
dev.off()

##################------------------------------
# Look for the variation in SLA in Room 1
png(file = "output/1.SLA_comparison.png", units="px", width=800, height=1200, res=200)
par(mfrow = c(3,1))
data.attrib.set = subset(data.attrib, Room %in% c(1)) # Consider the free seedling to test the parameter sensitivity

data.attrib.set.daily = summaryBy(Tgrowth+SLA+LA+Leafmass+Stemmass+Rootmass+g1+EaV+delsV+EaJ+delsJ+JVr+alpha+theta+Vcmax25+Jmax25+k+Y+af+as+ar+
                                    Intercept+Slope+self_s.mean ~ Room+Date, data=data.attrib.set, FUN=c(mean), na.rm=TRUE)
names(data.attrib.set.daily) = c("Room","Date","Tgrowth","SLA","LA","LM","SM","RM","g1","EaV","delsV","EaJ","delsJ","JVr","alpha","theta",
                                 "Vcmax25","Jmax25","k","Y","af","as","ar","Intercept","Slope","self_s")

plot(data.attrib.set.daily$Date,data.attrib.set.daily$SLA,col='red',type='l',main='Room 1',xlab='',ylab='SLA (m2 gC-1)')
lines(data.attrib.set.daily$Date,data.attrib.set.daily$LA/data.attrib.set.daily$LM,col='green')
legend('topleft', c("SLA from harvest", "Modelled SLA (allometry-based LA / modelled LM)"), lty=1, col=c('red','green'), bty='n', cex=0.75)

# Look for the variation in SLA in Room 4
data.attrib.set = subset(data.attrib, Room %in% c(4)) # Consider the free seedling to test the parameter sensitivity

data.attrib.set.daily = summaryBy(Tgrowth+SLA+LA+Leafmass+Stemmass+Rootmass+g1+EaV+delsV+EaJ+delsJ+JVr+alpha+theta+Vcmax25+Jmax25+k+Y+af+as+ar+
                                    Intercept+Slope+self_s.mean ~ Room+Date, data=data.attrib.set, FUN=c(mean), na.rm=TRUE)
names(data.attrib.set.daily) = c("Room","Date","Tgrowth","SLA","LA","LM","SM","RM","g1","EaV","delsV","EaJ","delsJ","JVr","alpha","theta",
                                 "Vcmax25","Jmax25","k","Y","af","as","ar","Intercept","Slope","self_s")

plot(data.attrib.set.daily$Date,data.attrib.set.daily$SLA,col='red',type='l',main='Room 4',xlab='',ylab='SLA (m2 gC-1)',ylim=c(0.04,0.06))
lines(data.attrib.set.daily$Date,data.attrib.set.daily$LA/data.attrib.set.daily$LM,col='green')
# legend('topleft', c("SLA from harvest", "Modelled SLA from allometry (LA/LM)"), lty=1, col=c('red','green'), bty='n', cex=0.75)

# Look for the variation in SLA in Room 6
data.attrib.set = subset(data.attrib, Room %in% c(6)) # Consider the free seedling to test the parameter sensitivity

data.attrib.set.daily = summaryBy(Tgrowth+SLA+LA+Leafmass+Stemmass+Rootmass+g1+EaV+delsV+EaJ+delsJ+JVr+alpha+theta+Vcmax25+Jmax25+k+Y+af+as+ar+
                                    Intercept+Slope+self_s.mean ~ Room+Date, data=data.attrib.set, FUN=c(mean), na.rm=TRUE)
names(data.attrib.set.daily) = c("Room","Date","Tgrowth","SLA","LA","LM","SM","RM","g1","EaV","delsV","EaJ","delsJ","JVr","alpha","theta",
                                 "Vcmax25","Jmax25","k","Y","af","as","ar","Intercept","Slope","self_s")

plot(data.attrib.set.daily$Date,data.attrib.set.daily$SLA,col='red',type='l',main='Room 6',xlab='',ylab='SLA (m2 gC-1)',ylim=c(0.04,0.08))
lines(data.attrib.set.daily$Date,data.attrib.set.daily$LA/data.attrib.set.daily$LM,col='green')
# legend('topleft', c("SLA from harvest", "Modelled SLA from allometry (LA/LM)"), lty=1, col=c('red','green'), bty='n', cex=0.75)

dev.off()

##################------------------------------
# take only Room 1 and 4 for attribution analysis
data.attrib = subset(data.attrib, Room %in% c(1,4))
met.rd.attrib = subset(met.rd.attrib, Room %in% c(1,4))

##################------------------------------
# Consider everything (Cday, LA, Rd, Biomass, parameters) for Room 1
q=0 # Case 0
data.attrib.set = subset(data.attrib, Room %in% c(1)) # Start with Room 4 to test the parameter sensitivity
# data.attrib.set$SLA = (subset(data.attrib, Room %in% c(4)))$SLA # Consider SLA from Room 6 to match the attribution
# data.attrib.set$LA = (subset(data.attrib, Room %in% c(4)))$LA # Consider LA from Room 6 to match the attribution
# data.attrib.set$Leafmass = (subset(data.attrib, Room %in% c(4)))$Leafmass # Consider LM from Room 6 to match the attribution

# data.attrib.test = subset(data.attrib.set, Date %in% as.Date("2016-01-08"))
# sum(data.attrib.test$R_leaf)

output.attrib.set = subset(output.attrib, Room %in% c(1))
met.rd.attrib.set = subset(met.rd.attrib, Room %in% c(1))

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# Meteorological data attribution
q=1 # Case 1
# Modify the Vcmax and Jmax for Room 4
# data.attrib.set$VPD = (subset(data.attrib, Room %in% c(4)))$VPD
data.attrib.set$Tair = (subset(data.attrib, Room %in% c(4)))$Tair
# data.attrib.set$PAR = (subset(data.attrib, Room %in% c(4)))$PAR

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# Meteorological data attribution
q=2 # Case 2
# Modify the Vcmax and Jmax for Room 4
data.attrib.set$VPD = (subset(data.attrib, Room %in% c(4)))$VPD
# data.attrib.set$Tair = (subset(data.attrib, Room %in% c(4)))$Tair
# data.attrib.set$PAR = (subset(data.attrib, Room %in% c(4)))$PAR

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# Meteorological data attribution
q=3 # Case 3
# Modify the respiration due to direct effect of temperature for warmed treatment
met.rd.attrib.set$Tair = (subset(met.rd.attrib, Room %in% c(4)))$Tair

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# Rm (leaf, stem and root)
q=4 # Case 4
# Modify the R_leaf, R_stem and R_root for Room 4
met.rd.attrib.set$R25_leaf = (subset(met.rd.attrib, Room %in% c(4)))$R25_leaf
met.rd.attrib.set$R25_stem = (subset(met.rd.attrib, Room %in% c(4)))$R25_stem
met.rd.attrib.set$R25_root= (subset(met.rd.attrib, Room %in% c(4)))$R25_root

# data.attrib.set$R_leaf = (subset(data.attrib, Room %in% c(4)))$R_leaf
# data.attrib.set$R_stem = ((subset(data.attrib, Room %in% c(4)))$R_stem)
# data.attrib.set$R_root = ((subset(data.attrib, Room %in% c(4)))$R_root)

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# Vcmax and Jmax attribution
q=5 # Case 5
# Modify the Vcmax and Jmax for Room 4
data.attrib.set$Vcmax25 = (subset(data.attrib, Room %in% c(4)))$Vcmax25
data.attrib.set$Jmax25 = (subset(data.attrib, Room %in% c(4)))$Jmax25

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# Temperature response parameters of Vcmax and Jmax (Activation energy (Ea), entropy (∆S))
q=6 # Case 6
# Modify the Ea and ∆S for Room 4
data.attrib.set$EaV = (subset(data.attrib, Room %in% c(4)))$EaV
data.attrib.set$delsC = (subset(data.attrib, Room %in% c(4)))$delsC
data.attrib.set$EaJ = (subset(data.attrib, Room %in% c(4)))$EaJ
data.attrib.set$delsJ = (subset(data.attrib, Room %in% c(4)))$delsJ

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# g1
q=7 # Case 7
# Modify the g1 for Room 4
data.attrib.set$g1 = (subset(data.attrib, Room %in% c(4)))$g1

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# alpha and theta
q=8 # Case 8
# Modify the alpha for Room 4
data.attrib.set$alpha = (subset(data.attrib, Room %in% c(4)))$alpha
data.attrib.set$theta = (subset(data.attrib, Room %in% c(4)))$theta

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# SLA
q=9 # Case 9
# Modify the SLA for Room 4
data.attrib.set$SLA = (subset(data.attrib, Room %in% c(4)))$SLA # Consider SLA from Room 4 to match the attribution
data.attrib.set$LA = (subset(data.attrib, Room %in% c(4)))$LA # Consider LA from Room 4 to match the attribution
data.attrib.set$Leafmass = (subset(data.attrib, Room %in% c(4)))$Leafmass # Consider LM from Room 4 to match the attribution

# data.attrib.set$SLA = (subset(data.attrib, Room %in% c(1)))$SLA # Consider SLA from Room 4 to match the attribution
# data.attrib.set$LA = (subset(data.attrib, Room %in% c(1)))$LA # Consider LA from Room 4 to match the attribution
# data.attrib.set$Leafmass = (subset(data.attrib, Room %in% c(1)))$Leafmass # Consider LM from Room 4 to match the attribution

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# af, as, ar
q=10 # Case 10
# Modify the af, as and ar for Room 4
data.attrib.set$af = (subset(data.attrib, Room %in% c(4)))$af
data.attrib.set$as = (subset(data.attrib, Room %in% c(4)))$as
data.attrib.set$ar = (subset(data.attrib, Room %in% c(4)))$ar

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# Y
q=11 # Case 11
# Modify the af, as and ar for Room 4
data.attrib.set$Y = (subset(data.attrib, Room %in% c(4)))$Y

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

##################------------------------------
# k
q=12 # Case 12
# Modify the af, as and ar for Room 4
data.attrib.set$k = (subset(data.attrib, Room %in% c(4)))$k

# This sript runs the model equations for parameter shifting from potted seedling to free seedling
source("R/CBM_model_shift_great.R")

# ##################------------------------------
# # self shading factors
# q=12 # Case 12
# # Modify the Slope and Intercept for Room 4
# data.attrib.set$Slope = (subset(data.attrib, Room %in% c(4)))$Slope
# data.attrib.set$Intercept = (subset(data.attrib, Room %in% c(4)))$Intercept
# 
# # This sript runs the model equations for parameter shifting from potted seedling to free seedling
# source("R/CBM_model_shift_great.R")
# 

############ Summarize the C pools plots
# Load the custom analysis and plotting functions that do all of the actual work	
source("R/functions_great.R")	

plot.shift = list() 
font.size = 7
title = as.character(c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)"))
cbPalette1 = c("gray", "skyblue")

data.attrib$Date = as.Date(data.attrib$Date)
data.attrib.daily = summaryBy(VPD+Tair+PAR+Tgrowth+SLA+LA+Leafmass+Stemmass+Rootmass+g1+EaV+delsV+EaJ+delsJ+JVr+alpha+theta+Vcmax25+Jmax25+k+Y+af+as+ar+
                                Intercept+Slope+self_s.mean ~ Room+Date, data=data.attrib, FUN=c(mean), na.rm=TRUE)
data.attrib.daily.resp = summaryBy(R_leaf+R_stem+R_root ~ Room+Date, data=data.attrib, FUN=c(sum), na.rm=TRUE)
data.attrib.daily = merge(data.attrib.daily, data.attrib.daily.resp, by=c("Room", "Date"))
names(data.attrib.daily) = c("Room","Date","VPD","Tair","PAR","Tgrowth","SLA","LA","LM","SM","RM","g1","EaV","delsV","EaJ","delsJ","JVr","alpha","theta",
                             "Vcmax25","Jmax25","k","Y","af","as","ar","Intercept","Slope","self_s","R_leaf","R_stem","R_root")

######## Plot Meteorological data
plot.shift[[1]] = plot.VPD.v1(data.attrib.daily, 1)
plot.shift[[2]] = plot.Tair(data.attrib.daily, 2)
# plot.shift[[3]] = plot.PAR(data.attrib.daily, 3)

######## Plot both Vcmax25 and Jmax25
plot.shift[[3]] = plot.Vcmax.Jmax.v1(data.attrib.daily, 3)
# plot.shift[[3]] = plot.Jmax(data.attrib.daily, 5)

######## Plot Respiration rates
plot.shift[[4]] = plot.Rd(data.attrib.daily, 4)

######## Plot Activation energy (Ea) and entropy (∆S)
# get Vcmax temperature response
tleaf<-seq(15,45,1)
data.attrib.daily.rm1 = subset(data.attrib.daily, Room == 1) 
room1 = fit_parr_vs(Ts=tleaf,Ea=mean(data.attrib.daily.rm1$EaV),delS=mean(data.attrib.daily.rm1$delsV))
Vcmax.response.rm1 = data.frame(tleaf=tleaf, Room = 1, response = room1)
data.attrib.daily.rm4 = subset(data.attrib.daily, Room == 4) 
room4 = fit_parr_vs(Ts=tleaf,Ea=mean(data.attrib.daily.rm4$EaV),delS=mean(data.attrib.daily.rm4$delsV))
Vcmax.response.rm4 = data.frame(tleaf=tleaf, Room = 4, response = room4)
Vcmax.response = rbind(Vcmax.response.rm1,Vcmax.response.rm4)

# plot.shift[[5]] = plot.Ea(data.attrib.daily, 1)
plot.shift[[5]] = plot.Ea(Vcmax.response, 5)

# get Jmax temperature response
data.attrib.daily.rm1 = subset(data.attrib.daily, Room == 1) 
room1 = fit_parr_vs(Ts=tleaf,Ea=mean(data.attrib.daily.rm1$EaJ),delS=mean(data.attrib.daily.rm1$delsJ))
Jmax.response.rm1 = data.frame(tleaf=tleaf, Room = 1, response = room1)
data.attrib.daily.rm4 = subset(data.attrib.daily, Room == 4) 
room4 = fit_parr_vs(Ts=tleaf,Ea=mean(data.attrib.daily.rm4$EaJ),delS=mean(data.attrib.daily.rm4$delsJ))
Jmax.response.rm4 = data.frame(tleaf=tleaf, Room = 4, response = room4)
Jmax.response = rbind(Jmax.response.rm1,Jmax.response.rm4)

# plot.shift[[5]] = plot.dels(data.attrib.daily, 1)
plot.shift[[6]] = plot.dels(Jmax.response, 6)

######## Plot g1
plot.shift[[7]] = plot.g1(data.attrib.daily, 7)

######## Plot alpha and theta
PAR<-seq(-100,2000,10)
data.attrib.daily.rm1 = subset(data.attrib.daily, Room == 1) 
room1 = curve.nlslrc(data.attrib.daily.rm1$alpha, data.attrib.daily.rm1$theta, mean(data.attrib.daily.rm1$Jmax25))
light.response.rm1 = data.frame(PAR=PAR, Room = 1, response = room1$J)
data.attrib.daily.rm4 = subset(data.attrib.daily, Room == 4) 
room4 = curve.nlslrc(data.attrib.daily.rm4$alpha, data.attrib.daily.rm4$theta, mean(data.attrib.daily.rm4$Jmax25))
light.response.rm4 = data.frame(PAR=PAR, Room = 4, response = room4$J)
light.response = rbind(light.response.rm1,light.response.rm4)

plot.shift[[8]] = plot.alpha.theta(light.response, 8)

######## Plot SLA
sla.harvest.all = read.csv("processed_data/sla.harvest.all.csv")
sla.harvest.all = subset(sla.harvest.all, Room %in% as.factor(c("1","4")))
sla.harvest.all$Date = as.Date(sla.harvest.all$Date)

plot.shift[[9]] = plot.SLA.v1(data.attrib.daily, sla.harvest.all, 9)

# Plot individual modelled parameters ("k","Y","af","sf") against "volume"
plot.shift[[10]] = plot.allocation.fractions(data.attrib.daily, 10)

plot.shift[[11]] = plot.Y(data.attrib.daily, 11)

plot.shift[[12]] = plot.k(data.attrib.daily, 12)


png("output/Figure_parameter_shifting_rm1-4_params.png", units="px", width=3000, height=1500, res=250)
lay <- rbind(c(1,2,3),c(4,5,6),c(7,8,9),c(10,11,12))
grid.arrange(grobs = plot.shift, layout_matrix = lay)
dev.off()

# # This sript plots the C pools for various test cases with parameter shifted from potted seedling to free seedling
# source("R/generate_figures_param_shift.R")
font.size = 12
# cbPalette = c("gray", "orange", "skyblue", "green", "black", "yellow3", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "mediumorchid", "lightpink", "olivedrab1")
cbPalette = c("gray", "orange", "skyblue", "green", "black", "yellow3", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "mediumorchid", "olivedrab1","navy")

shift.output.Mleaf = subset(shift.output,(variable %in% "Mleaf"))
plot.shift[[13]] = plot.Mleaf.rm1to4(shift.output.Mleaf)

shift.output.Mstem = subset(shift.output,(variable %in% "Mstem"))
plot.shift[[14]] = plot.Mstem.rm1to4(shift.output.Mstem)

shift.output.Mroot = subset(shift.output,(variable %in% "Mroot"))
plot.shift[[15]] = plot.Mroot.rm1to4(shift.output.Mroot)

png("output/Figure_parameter_shifting_rm1-4_linegraphs.png", units="px", width=3000, height=3500, res=300)
lay <- rbind(c(1,2,13,13),c(3,4,13,13),c(5,6,14,14),c(7,8,14,14),c(9,10,15,15),c(11,12,15,15))
grid.arrange(grobs = plot.shift, layout_matrix = lay)
dev.off()

shift.output.dcast = reshape2::dcast(shift.output, Date+Case ~ variable)
shift.output.dcast$biomass = shift.output.dcast$Mleaf + shift.output.dcast$Mstem + shift.output.dcast$Mroot
# plot.shift[[16]] = plot.biomass.rm1to4(shift.output.dcast)
plot.shift[[16]] = plot.biomass.logscale.rm1to4(shift.output.dcast)

shift.output.bar = subset(shift.output.dcast, Date %in% as.Date("2016-02-24"))
keeps = c("Date","Case","biomass")
shift.output.bar = shift.output.bar[ , keeps, drop = FALSE]
# shift.output.bar[,c("Date","Case","Cstorage","Sleaf")] = NULL
row.names(shift.output.bar) <- as.numeric(paste(shift.output.bar$Case))
# shift.output.bar[,c("Date","Case")] = NULL
plot.shift[[17]] = plot.biomass.barplot.rm1to4(shift.output.bar)
plot1 = plot.shift[[17]] +
  annotate("text", x = max(as.numeric(shift.output.bar$Case)), y = max(shift.output.bar$biomass), size = font.size-7, label = paste("(a)"))

png("output/Figure_parameter_shifting_rm1-4_barchart.png", units="px", width=1000, height=800, res=200)
plot.shift[[17]]
dev.off()

png("output/Figure_parameter_shifting_rm1-4.png", units="px", width=3000, height=4000, res=300)
lay <- rbind(c(1,2,13,13),c(3,4,13,13),c(5,6,14,14),c(7,8,14,14),c(9,10,15,15),c(11,12,15,15),c(16,16,17,17),c(16,16,17,17))
grid.arrange(grobs = plot.shift, layout_matrix = lay)
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


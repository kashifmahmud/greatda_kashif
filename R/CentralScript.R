#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- Central analysis script for DA with GREAT experiment manuscript, entitled
#  "".

#  The idea is to keep this script nice and tidy, but reproducibly do all the
#  analysis and make all of the figures for the manuscript. Raw and processed data will be 
#  placed in the "data" and "processed_data" folders respectively, while figures 
#  and tables will be placed in the "output" folder.
#-------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------
# Developed by Kashif Mahmud (March 2018)
# k.mahmud@westernsydney.edu.au
#----------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# # Clear the workspace (if needed)
# rm(list=ls())
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# Load required libraries. There are quite a few, including some non-standard functions that are not on CRAN. 
# This script will check for required libraries and install any that are missing	
source('R/load_packages_great.R')	

# Load the custom analysis and plotting functions that do all of the actual work	
source("R/functions_great.R")	
source("R/functions_great_CBM.R")	
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# #- Download data files for WTC3 experiment. This downloads the zipfile from figshare
# download.file("https://ndownloader.figshare.com/files/4857112", "data.zip", mode="wb")
# # Extract data to different folders.
# unzip("data.zip")

#-------------------------------------------------------------------------------------

# #-------------------------------------------------------------------------------------
#- This script imports and processes the raw GREAT experiment biomass data to model the carbon pools and fluxes using DA
# source("R/initial_biomass_data_processing.temperature.R")
# Script modification to add harvest Height and Diameter to estimate the harvest biomass to end the experiment on harvest
source("R/initial_biomass_data_processing.temperature_v2.R")

# rmd2rscript("report_initial_data_processing_wtc3.Rmd")
# source("report_initial_data_processing_wtc3.R")
# #-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Make figure 1. Model representation of storage, allocation and autotrophic respiration processes and
# pathways in the CBM with storage pool, separate growth and maintenance respiration components.
plot.model.great()
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Read the processed data and clean the workspace from data pre-processing
#  Processed data are placed in "processed_data" folder
# source("R/read_data_great.R")

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# Model run using GREAT dataset (investigate the temperature effects by considering only the well-watered treatments)
# Read the pre-processed data of GPP, Respiration and biomass
data.biomass = read.csv("processed_data/modelled_data.csv") # Estimated biomass from height and dia
names(data.biomass) = c("Date","Room","LA","LA_SE","LM","LM_SE","WM","WM_SE","RM","RM_SE")
data.biomass$Date = as.Date(data.biomass$Date, format = "%Y-%m-%d")

# # Data up to the end of experiment (29th February)
# data.biomass = subset(data.biomass,(Date %in% as.Date(c("2016-01-08","2016-01-18","2016-01-28","2016-02-08","2016-02-29"))))
# data.gpp = read.csv("processed_data/great_daily_carbon_gain_LA_29-02-2016.csv")

# Data up to the final harvest (23rd February)
data.biomass = subset(data.biomass,(Date %in% as.Date(c("2016-01-08","2016-01-18","2016-01-28","2016-02-08","2016-02-23"))))
# data.gpp = read.csv("processed_data/great_daily_carbon_gain_LA_23-02-2016.csv")
data.gpp = read.csv("processed_data/great_daily_carbon_gain_LA_2016-02-24.csv") # Estimated biomass from height and dia
data.gpp$Date = as.Date(data.gpp$Date, format="%d/%m/%Y")

keeps = c("Date","Room","GPP","R_leaf","R_leaf_se","R_stem","R_stem_se","R_root","R_root_se")
data.gpp = data.gpp[ , keeps, drop = FALSE]
names(data.gpp) = c("Date","Room","GPP","R_leaf","R_leaf_SE","R_wood","R_wood_SE","R_root","R_root_SE")
data.gpp$Date = as.Date(data.gpp$Date, format = "%Y-%m-%d")
data.all = merge(data.gpp, data.biomass, by=c("Date","Room"), all=TRUE)

#-------------------------------------------------------------------------------------
# # Data up to the final harvest (23rd February)
# # data.biomass = read.csv("processed_data/harvest_data.csv") # Direct harvest data
# data.biomass = read.csv("processed_data/modelled_data.csv") # Estimated biomass data
# names(data.biomass) = c("Date","Room","LA","LA_SE","LM","LM_SE","WM","WM_SE","RM","RM_SE")
# data.biomass$Date = as.Date(data.biomass$Date, format = "%Y-%m-%d")
# data.biomass = data.biomass[data.biomass$Date < as.Date("2016-02-29"),]
# 
# data.gpp = read.csv("processed_data/great_daily_carbon_gain_LA_23-02-2016.csv")
# keeps = c("Date","Room","GPP","R_leaf","R_leaf_se","R_stem","R_stem_se","R_root","R_root_se")
# data.gpp = data.gpp[ , keeps, drop = FALSE]
# names(data.gpp) = c("Date","Room","GPP","R_leaf","R_leaf_SE","R_wood","R_wood_SE","R_root","R_root_SE")
# data.gpp$Date = as.Date(data.gpp$Date, format = "%Y-%m-%d")
# data.all = merge(data.gpp, data.biomass, by=c("Date","Room"), all=TRUE)
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# Plot Respiration rates
cbPalette = rainbow(6)[rank(1:6)]
p0 = ggplot(data.all, aes(x=Date, y=R_leaf, group = Room, colour=as.factor(Room))) + 
  geom_point(size=1) +
  geom_line(data = data.all, aes(x = Date, y = R_leaf, group = Room, colour=as.factor(Room))) +
  ylab(expression(Foliage~respiration~(g~C~g~C^{-1}~d^{-1}))) + 
  # xlab("Month") +
  # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g DM)")) +
  # ggtitle("C pools - Measured (points) vs Modelled (lines)") +
  labs(colour="Room") + scale_color_manual(values=cbPalette) +
  theme_bw() +
  # annotate("text", x = min(summary.error.Cpool$Date), y = max(summary.error.Cpool$value), size = 14, label = paste(title[p])) +
  # theme(plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.title = element_text(colour="black", size=12)) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  # theme(legend.key.height=unit(0.9,"line")) +
  theme(legend.position = c(0.2,0.7)) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 14, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

p0 = ggplot(data.all, aes(x=Date, y=R_wood, group = Room, colour=as.factor(Room))) + 
  geom_point(size=1) +
  geom_line(data = data.all, aes(x = Date, y = R_wood, group = Room, colour=as.factor(Room))) +
  ylab(expression(Wood~respiration~(g~C~g~C^{-1}~d^{-1}))) + 
  # xlab("Month") +
  # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g DM)")) +
  # ggtitle("C pools - Measured (points) vs Modelled (lines)") +
  labs(colour="Room") + scale_color_manual(values=cbPalette) +
  theme_bw() +
  # annotate("text", x = min(summary.error.Cpool$Date), y = max(summary.error.Cpool$value), size = 14, label = paste(title[p])) +
  # theme(plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.title = element_text(colour="black", size=12)) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  # theme(legend.key.height=unit(0.9,"line")) +
  theme(legend.position = c(0.2,0.7)) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 14, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

p0 = ggplot(data.all, aes(x=Date, y=R_root, group = Room, colour=as.factor(Room))) + 
  geom_point(size=1) +
  geom_line(data = data.all, aes(x = Date, y = R_root, group = Room, colour=as.factor(Room))) +
  ylab(expression(Root~respiration~(g~C~g~C^{-1}~d^{-1}))) + 
  # xlab("Month") +
  # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g DM)")) +
  # ggtitle("C pools - Measured (points) vs Modelled (lines)") +
  labs(colour="Room") + scale_color_manual(values=cbPalette) +
  theme_bw() +
  # annotate("text", x = min(summary.error.Cpool$Date), y = max(summary.error.Cpool$value), size = 14, label = paste(title[p])) +
  # theme(plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.title = element_text(colour="black", size=12)) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  # theme(legend.key.height=unit(0.9,"line")) +
  theme(legend.position = c(0.2,0.7)) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 14, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

pdf(file = "output/Respiration_rates_DA.pdf",width=10, height=5)
print (p0)
dev.off() 
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# Plot GPP
cbPalette = rainbow(6)[rank(1:6)]
p0 = ggplot(data.all, aes(x=Date, y=GPP, group = Room, colour=as.factor(Room))) + 
  geom_point(size=1) +
  geom_line(data = data.all, aes(x = Date, y = GPP, group = Room, colour=as.factor(Room))) +
  ylab(expression(GPP~"("*g~C~d^"-1"*")")) + 
  # xlab("Month") +
  # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g DM)")) +
  # ggtitle("C pools - Measured (points) vs Modelled (lines)") +
  labs(colour="Room") + scale_color_manual(values=cbPalette) +
  theme_bw() +
  # annotate("text", x = min(summary.error.Cpool$Date), y = max(summary.error.Cpool$value), size = 14, label = paste(title[p])) +
  # theme(plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.title = element_text(colour="black", size=12)) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  # theme(legend.key.height=unit(0.9,"line")) +
  theme(legend.position = c(0.2,0.7)) +
  theme(legend.key = element_blank()) +
  theme(text = element_text(size=12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 14, vjust=0.3)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# png(file = "Output/GPP.png")
pdf(file = "output/GPP_v2.pdf",width=10, height=5)
print (p0)
dev.off() 
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Matching C balance of the entire experiment considering C inputs and outputs
source("R/C_balance_great.R")

#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# Calculate SLA from data
source("R/sla_great.R")

#-------------------------------------------------------------------------------------
# Load required libraries and functions in case you clear the workspace after pre-processing
source("R/functions_great.R")
source("R/functions_great_CBM.R")

# Assign treatment groups
treat.group = unique(as.factor(data.all$Room)) # Assign all treatments
# treat.group = as.factor(c("1","4","6")) # Assign few treatments to check the results

# Data up to the 24th February
data.all = subset(data.all, Date <= as.Date("2016-02-24"))

# Model run for WTC3 dataset with clustering
cluster <- makeCluster(detectCores()-1)
# clusterEvalQ(cluster, library(xts))
clusterExport(cl=cluster, list("data.all","treat.group","tnc"))
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterExport(cluster, ex)
result.cluster = list()
bic.cluster = list()

start <- proc.time() # Start clock
# result <- clusterMap(cluster, mcmc.great, with.storage=rep(T,6), model.comparison=rep(F,6), model.optimization=rep(F,6), 
#                      no.param.par.var=rep(3,6),
#                      treat.group=treat.group,
#                      MoreArgs=list(chainLength=3000))
result <- clusterMap(cluster, mcmc.great, treat.group=treat.group,
                     MoreArgs=list(chainLength=11111,with.storage=T, model.comparison=F, 
                     model.optimization=F, no.param.per.var=2))

time_elapsed_series <- proc.time() - start # End clock
stopCluster(cluster)

listOfDataFrames <- vector(mode = "list", length = nlevels(treat.group))
for (i in 1:nlevels(treat.group)) {
  listOfDataFrames[[i]] <- data.frame(result[[i]][[6]])
}
bic = do.call("rbind", listOfDataFrames)
write.csv(bic, "output/bic.csv", row.names=FALSE)
sum(bic$bic)

# Plot parameters and biomass data fit
plot.Modelled.parameters.great(result,with.storage=T,treat.group)
plot.Modelled.biomass.great(result,with.storage=T,treat.group)

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# # Load required libraries and functions in case you clear the workspace after pre-processing
# source("R/functions_great.R")	
# source("R/functions_great_CBM.R")	
# start <- proc.time() # Start clock
# # 3000 chain length is sufficient for the convergance
# 
# # Provide inputs for the MCMC simulation
# chainLength = 2000 # Try larger iteration (chain length) if the mcmc chain doesn't converge (Chech the loglikelihood convergence with chain length)
# no.param.par.var = 2 # Parameter setting: Linear=2 / Quadratic=3 / Cubic=4
# with.storage = T # (default with storage)
# model.comparison = F # (default)
# model.optimization = F # (default)
# 
# # Main function to run MCMC simulation
# result = mcmc.great(chainLength, no.param.par.var, treat.group, with.storage, model.comparison, model.optimization) # Linear/Quadratic/Cubic parameters
# # result = CBM.wtc3(chainLength = 3000, no.param.par.var=(nrow(data.all)/4)/30, treat.group=treat.group, with.storage, model.comparison=F, model.optimization=F) # Monthly parameters
# 
# time_elapsed_series <- proc.time() - start # End clock
# 
# result[[6]]
# write.csv(result[[6]], "output/bic.csv", row.names=FALSE) # unit of respiration rates: gC per gC plant per day	

#-------------------------------------------------------------------------------------
# # Plot parameters and biomass data fit
# plot.Modelled.parameters(result,with.storage=T)
# plot.Modelled.biomass(result,with.storage=T)
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
# Calculate total C partitioning for individual treatments 
# and make figure 7 and Table S1
source("R/C_partitioning_great.R")
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Perform sensitivity analysis and make figure 5
# source("R/functions_CBM.R")
source("R/Param_sensitivity_great_rm1-4.R")
source("R/Param_sensitivity_great_rm4-6.R")

source("R/Parameter_sensitivity_great_rm1-4_individual_param.R")
source("R/Parameter_sensitivity_great_rm4-1_individual_param.R")
source("R/Parameter_sensitivity_great_rm4-6_individual_param.R")
source("R/Parameter_sensitivity_great_rm6-4_individual_param.R")
#-------------------------------------------------------------------------------------

# #-------------------------------------------------------------------------------------
# 
# # Check the C input and output balance
# data.set = subset(data.all,(Treatment %in% treat.group[v]))
# input = sum(data.set$GPP)
# Rm = subset(result[[4]], variable %in% as.factor("Rm"))
# Rm.sum = sum(Rm$value)
# 
# Rm.daily = mean(data.set$Rd.foliage.mean)*(data.set$LM[1]+data.set$LM[nrow(data.set)])/2 + mean(data.set$Rd.stem.mean)*(data.set$WM[1]+data.set$WM[nrow(data.set)])/2*mean(data.set$SMratio) + 
#   mean(data.set$Rd.branch.mean)*(data.set$WM[1]+data.set$WM[nrow(data.set)])/2*mean(data.set$BMratio) + 
#   mean(data.set$Rd.fineroot.mean)*(data.set$RM[1]+data.set$RM[nrow(data.set)])/2*mean(data.set$FRratio_SE) + mean(data.set$Rd.intermediateroot.mean)*(data.set$RM[1]+data.set$RM[nrow(data.set)])/2*mean(data.set$IRratio) + 
#   mean(data.set$Rd.coarseroot.mean)*(data.set$RM[1]+data.set$RM[nrow(data.set)])/2*mean(data.set$CRratio_SE) + mean(data.set$Rd.boleroot.mean)*(data.set$RM[1]+data.set$RM[nrow(data.set)])/2*mean(data.set$BRratio) 
# Rm.sum = Rm.daily * 252
#   
# output = Rm.sum + (data.set$LM[nrow(data.set)] - data.set$LM[1]) + (data.set$WM[nrow(data.set)] - data.set$WM[1]) + (data.set$RM[nrow(data.set)] - data.set$RM[1]) + 
#   (data.set$TNC_tot[max(which(complete.cases(data.set$TNC_tot)))] - data.set$TNC_tot[min(which(complete.cases(data.set$TNC_tot)))]) + data.set$litter[max(which(complete.cases(data.set$litter)))]
# 
# Y.modelled = subset(result[[2]], variable %in% as.factor("Y"))
# 
# #-------------------------------------------------------------------------------------
# # Check Rabove from data
# Rabove.ini = data.set$Rd.foliage.mean[1]*data.set$LM[1] + data.set$Rd.stem.mean[1]*data.set$WM[1]*data.set$SMratio[1] + data.set$Rd.branch.mean[1]*data.set$WM[1]*data.set$BMratio[1] +
#   Y.modelled$Parameter[1]*((data.set$LM[nrow(data.set)]-data.set$LM[1])/251 + (data.set$WM[nrow(data.set)]-data.set$WM[1])/251)
# 
# Rabove.end = data.set$Rd.foliage.mean[nrow(data.set)]*data.set$LM[nrow(data.set)] + data.set$Rd.stem.mean[nrow(data.set)]*data.set$WM[nrow(data.set)]*data.set$SMratio[nrow(data.set)] + data.set$Rd.branch.mean[nrow(data.set)]*data.set$WM[nrow(data.set)]*data.set$BMratio[nrow(data.set)] +
#     Y.modelled$Parameter[nrow(data.set)]*((data.set$LM[nrow(data.set)]-data.set$LM[1])/251 + (data.set$WM[nrow(data.set)]-data.set$WM[1])/251)
# 
# # mean(data.set$Ra)/mean(data.set$GPP)
# 
# # Rabove = mean(data.set$Rd.foliage.mean)*(data.set$LM[1]+data.set$LM[nrow(data.set)])/2 + mean(data.set$Rd.stem.mean)*(data.set$WM[1]+data.set$WM[nrow(data.set)])/2*mean(data.set$SMratio)
# # + mean(data.set$Rd.branch.mean)*(data.set$WM[1]+data.set$WM[nrow(data.set)])/2*mean(data.set$BMratio) +
# #   mean(Y.modelled$Parameter)*(data.set$LM[nrow(data.set)]-data.set$LM[1] + data.set$WM[nrow(data.set)]-data.set$WM[1])
# # Ra = sum(data.set$Ra)
# 
# #-------------------------------------------------------------------------------------
# 
# #-------------------------------------------------------------------------------------
# 
# cluster <- makeCluster(detectCores()-1)
# # clusterEvalQ(cluster, library(xts))
# clusterExport(cl=cluster, list("treat.group","data.all","tnc.partitioning"))
# ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
# clusterExport(cluster, ex)
# result.cluster = list()
# bic.cluster = list()
# 
# start <- proc.time() # Start clock
# result.cluster <- clusterMap(cluster, CBM.wtc3, no.param.par.var=c((nrow(data.all)/4)/30,(nrow(data.all)/4)/7), 
#             MoreArgs=list(chainLength=100, treat.group=treat.group, with.storage=c(T,T,T,T), model.comparison=c(F,F,F,F), model.optimization=c(F,F,F,F)))
# 
# time_elapsed_series <- proc.time() - start # End clock
# bic.without.storage = result.cluster[[1]][[6]]
# bic.with.storage = result.cluster[[2]][[6]]
# bic.group1 = result.cluster[[3]]
# bic.group2 = result.cluster[[4]]
# bic.group3 = result.cluster[[5]]
# bic.group4 = result.cluster[[6]]
# result = result.cluster[[7]]
# stopCluster(cluster)
# 
# #-------------------------------------------------------------------------------------
# 
# #-------------------------------------------------------------------------------------




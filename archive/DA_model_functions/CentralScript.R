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
rm(list=ls())
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# Load required libraries. There are quite a few, including some non-standard functions that are not on CRAN. 
# This script will check for required libraries and install any that are missing	
source('R/DA_model_functions/load_packages_great.R')	

# Load the custom analysis and plotting functions that do all of the actual work	
source("R/DA_model_functions/functions_great.R")	
source("R/DA_model_functions/functions_great_CBM.R")	
#-------------------------------------------------------------------------------------


# #-------------------------------------------------------------------------------------
#- This script imports and processes the raw GREAT experiment biomass data to model the carbon pools and fluxes using DA
# source("R/initial_biomass_data_processing.temperature.R")
# Script modification to add harvest Height and Diameter to estimate the harvest biomass to end the experiment on harvest

source("R/DA_model_functions/initial_biomass_data_processing.temperature_v2.R")

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

# Model run using GREAT dataset (investigate the temperature effects by considering only the well-watered treatments)
# Read the pre-processed data of GPP, Respiration and biomass

data.biomass = read.csv("processed_data/modelled_data.csv") # Estimated biomass from height and dia
names(data.biomass) = c("Date","Room","LA","LA_SE","LM","LM_SE","WM","WM_SE","RM","RM_SE")
data.biomass$Date = as.Date(data.biomass$Date, format = "%Y-%m-%d")

# Data up to the end of experiment (29th February)
# data.biomass = subset(data.biomass,(Date %in% as.Date(c("2016-01-08","2016-01-18","2016-01-28","2016-02-08","2016-02-29"))))
# data.gpp = read.csv("processed_data/great_daily_carbon_gain_LA_2016-02-29.csv")

# # Data up to the final harvest (24th February)

data.biomass = subset(data.biomass,(Date %in% as.Date(c("2016-01-08","2016-01-18","2016-01-28","2016-02-08","2016-02-23"))))
data.gpp = read.csv("processed_data/great_daily_carbon_gain_LA_2016-02-24.csv")

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
# # Plot GPP
# cbPalette = rainbow(6)[rank(1:6)]
# p0 = ggplot(data.all, aes(x=Date, y=GPP, group = Room, colour=as.factor(Room))) + 
#   geom_point(size=1) +
#   geom_line(data = data.all, aes(x = Date, y = GPP, group = Room, colour=as.factor(Room))) +
#   ylab(expression(GPP~"("*g~C~d^"-1"*")")) + 
#   # xlab("Month") +
#   # coord_trans(y = "log10") + ylab(paste(as.character(meas[p]),"(g DM)")) +
#   # ggtitle("C pools - Measured (points) vs Modelled (lines)") +
#   labs(colour="Room") + scale_color_manual(values=cbPalette) +
#   theme_bw() +
#   # annotate("text", x = min(summary.error.Cpool$Date), y = max(summary.error.Cpool$value), size = 14, label = paste(title[p])) +
#   # theme(plot.title = element_text(size = 20, face = "bold")) +
#   theme(legend.title = element_text(colour="black", size=12)) +
#   theme(legend.text = element_text(colour="black", size = 12)) +
#   # theme(legend.key.height=unit(0.9,"line")) +
#   theme(legend.position = c(0.2,0.7)) +
#   theme(legend.key = element_blank()) +
#   theme(text = element_text(size=12)) +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.title.y = element_text(size = 14, vjust=0.3)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
# 
# # png(file = "Output/GPP.png")
# pdf(file = "Output/GPP_v2.pdf",width=10, height=5)
# print (p0)
# dev.off() 
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- Matching C balance of the entire experiment considering C inputs and outputs
source("R/DA_model_functions/C_balance_great.R")

#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# Calculate SLA from data
source("R/DA_model_functions/sla_great.R")

#-------------------------------------------------------------------------------------
# Load required libraries and functions in case you clear the workspace after pre-processing
source("R/DA_model_functions/functions_great.R")
source("R/DA_model_functions/functions_great_CBM.R")

# Assign treatment groups
 treat.group = unique(as.factor(data.all$Room)) # Assign all treatments
# treat.group = as.factor(c("1","4","6")) # Assign few treatments to check the results

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

result <- clusterMap(cluster, mcmc.great, 
                     treat.group=treat.group,
                     MoreArgs=list(chainLength=5000,with.storage=T, model.comparison=F, model.optimization=F, 
                                   no.param.per.var=2))

time_elapsed_series <- proc.time() - start # End clock
stopCluster(cluster)

listOfDataFrames <- vector(mode = "list", length = nlevels(treat.group))
for (i in 1:nlevels(treat.group)) {
  listOfDataFrames[[i]] <- data.frame(result[[i]][[6]])
}


bic = do.call("rbind", listOfDataFrames)
write.csv(bic, "output/bic.csv", row.names=FALSE)
sum(bic$bic)


#-------------------------------------------------------------------------------------

# Plot parameters and biomass data fit

source("R/plot_DA_results.R")

plot_da_parameters(result=result)
plot_da_nsc(result=result)
plot_da_data(result=result)
plot_da_data_log(result=result)
#-------------------------------------------------------------------------------------

# plot.Modelled.parameters.great(result,with.storage=T,treat.group)
# plot.Modelled.biomass.great(result,with.storage=T,treat.group)
# 


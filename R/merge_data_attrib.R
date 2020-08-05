# merge attribution data (from Dushan) and new DA results

listOfDataFrames <- vector(mode = "list", length = nlevels(treat.group))
for (i in 1:nlevels(treat.group)) {
  listOfDataFrames[[i]] <- data.frame(result[[i]][[2]])
}
summary.param = do.call("rbind", listOfDataFrames)

keeps = c("Date", "variable", "Parameter", "treatment")
param.attrib = summary.param[ , keeps, drop = FALSE]
param.attrib = reshape2::dcast( param.attrib , Date+treatment ~ variable, value.var = "Parameter")
names(param.attrib)[2] = "Room"


listOfDataFrames <- vector(mode = "list", length = nlevels(treat.group))
for (i in 1:nlevels(treat.group)) {
  listOfDataFrames[[i]] <- data.frame(result[[i]][[4]])
}
summary.output = do.call("rbind", listOfDataFrames)

keeps = c("Date", "variable", "value", "treatment")
output.attrib = summary.output[ , keeps, drop = FALSE]
output.attrib = reshape2::dcast( output.attrib , Date+treatment ~ variable, value.var = "value")
names(output.attrib) = c("Date","Room",'Leafmass','Stemmass','Rootmass',"Rm")

#--------------------------------------
# Read tree attributes data including Met data, Temperature dependant variables, Modelled Parameters from Dushan's analysis 
# data.attrib <- read.csv("parameters/data_for_attribution_analysis_v2.csv")
# data.attrib$DateTime_hr = as.POSIXct(as.character(data.attrib$DateTime_hr), format="%d/%m/%Y %H:%M")
# data.attrib$Date = as.Date(data.attrib$Date, format="%d/%m/%Y")
data.attrib <- read.csv("parameters/data_for_attribution_analysis_v3.csv")
data.attrib$Tgrowth[data.attrib$Tgrowth == 25.8] = 28.5 #correction to temperature from 25.8 to 28.5 
data.attrib$DateTime_hr = as.POSIXct(as.character(data.attrib$DateTime_hr), format="%Y-%m-%d %H:%M")
data.attrib$Date = as.Date(data.attrib$Date, format="%Y-%m-%d")
data.attrib$Room = as.factor(data.attrib$Room)

# replace modelled biomass and DA parameters
# data.attrib[ ,c('k','Y','af','as','ar','Leafmass','Stemmass','Rootmass','Rm')] <- list(NULL)
data.attrib[ ,c('Leafmass','Stemmass','Rootmass')] <- list(NULL)
data.attrib = merge(data.attrib,param.attrib,by=c("Date","Room"))
data.attrib = merge(data.attrib,output.attrib,by=c("Date","Room"))

write.csv(data.attrib, "Parameters/data.attrib.merged.csv", row.names=FALSE)

#--------------------------------------





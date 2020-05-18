############################################################################################################
#one point Vcmax estimation
############################################################################################################

#some functions

# function to get Rgas
.Rgas <- function()8.314

# function to convert C to K
Tk <- function(x)x+273.15

# function to get Arrhenius model

arrh <- function(Tleaf, Ea){
  exp((Ea * (Tk(Tleaf) - 298.15)) / (298.15 * .Rgas() * Tk(Tleaf))) 
}


# function to get temperature response of gammastar (Bernachchi et al 2001 parameters)

TGammaStar <- function(Tleaf, Patm=100,
                       Egamma=37830.0, 
                       value25=42.75){  
  
  value25*arrh(Tleaf,Egamma)*Patm/100
}


#function to get temperature response of Km (Bernachchi et al 2001 parameters)

T_Km <- function(Tleaf, Patm=100,
                 Oi = 210,      # O2 concentration (mmol mol-1)
                 Ec=79403,  # activation energy for Kc 
                 Eo=36380,  # activation energy for Ko
                 Kc25=404.9,  # Kc at 25C
                 Ko25=278.4 # Ko at 25C
){
  
  Oi <- Oi * Patm / 100
  
  Ko <- Ko25*arrh(Tleaf, Eo)
  Kc <- Kc25*arrh(Tleaf, Ec)
  Km <- Kc * (1.0 + Oi / Ko)
  
  return(Km)
}

# test functions
T_Km(30)
TGammaStar(25)

# function to get one point Vcmax

o_p_vcmax<-function(Tleaf,Ci,Photo){
  
  gstar<-TGammaStar(Tleaf)
  km<-T_Km(Tleaf)
  a1<-Ci-gstar
  a2<-Ci+km
  
  vcmax<-Photo/((a1/a2)-0.015)
  
  return(vcmax)
}

# o_p_vcmax(25,400,25)

data<-getAvT()
data$Vcmax<-o_p_vcmax(data$Tleaf,data$Ci,data$Photo)
with(data,plot(Tleaf,Vcmax))
summaryBy(Vcmax~Room,data=subset(data,data$PARi>1200),FUN=c(mean,std.error))

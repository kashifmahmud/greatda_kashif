

# function to fit peaked arrhenius model
fit_parr_vs<-function(Ts,Ea,delS){
  k25=1
  Vcmax<-k25 * exp((Ea*(Tk(Ts) - 298.15))/(298.15*0.008314*Tk(Ts))) * 
    (1+exp((298.15*delS - 200)/(298.15*0.008314))) / 
    (1+exp((Tk(Ts)*delS-200)/(Tk(Ts)*0.008314)))
  
  return(Vcmax)
}


# get Vcmax temperature response

tleaf<-seq(15,45,1)

room1_dat<-subset(data,data$Room==1) # subset of data for room 1....
room1<-fit_parr_vs(Ts=tleaf,Ea=mean(room1_dat$EaV),delS=mean(room1_dat$delsV))



# get Jmax temperature response
room1<-fit_parr_vs(Ts=tleaf,Ea=mean(room1_dat$EaJ),delS=mean(room1_dat$delsJ))

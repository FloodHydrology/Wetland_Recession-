#Greensboro------------------------------------------------------------------
gage<-'01491000'
df<-q_gage %>% filter(site_no == gage) %>%
  select(Date, X_00060_00003) %>%
  rename(Q = X_00060_00003) %>%
  na.omit()

#Convert to mm/day
df$Q<-df$Q/(gages$Watershed_Area[gages$SOURCE_FEA==gage])*(0.3048^3)/10*86400


if(nrow(df)>0){
  #Estimate baseflow using EcoHydRology package
  df$Q_bf <- BaseflowSeparation(df$Q, filter_parameter = 0.95)$bt
  
  #Estimate dQ/dt
  df<- df %>%
    mutate(dQ = Q_bf - lag(Q_bf)) %>%
    filter(dQ<0) %>%
    filter(Q_bf>0)
  
  # Take monthly means
  #df %>% group_by(month(Date)) %>% summarise(dQ=mean(dQ), Q=mean(Q)) %>% plot()  
  
  #Create model to esitmate a and b
  if(nrow(df)>0){
    model<-lm(log10(dQ*-1)~log10(Q), data=df)
    b<-model$coefficients[1]
    a<-model$coefficients[2]
    p_value<-summary(model)$coefficients[8]
  }else{
    a<- -9999
    b<- -9999
    p_value<- -9999
  }
}else{
  a<- -9999
  b<- -9999
  p_value<- -9999
}
plot(df$Q, df$dQ*-1, log="xy", 
     ylim=c(0.000001, 0.1), xlim=c(0.01, 10),
     ylab="dq/dt [mm/day]", xlab="Streamflow (mm)",
     cex.lab=14/12, cex.axis=10/12, 
     pch=19, col="#D3D3D366", main="Greensboro")
plot(df$Date, df$Q[617:649],df$dQ[617:649]*-1, type="l", col="dark red", lwd=2)
points(df$Q[1055:1090], df$dQ[1055:1090]*-1, type="l", lwd=2, col="dark blue")
points(df$Q[2395:2400], df$dQ[2395:2400]*-1, type="l", lwd=2, col="orange")

#Tuckahoe--------------------------------------------------------------------------
gage<-'01491500'
df<-q_gage %>% filter(site_no == gage) %>%
  select(Date, X_00060_00003) %>%
  rename(Q = X_00060_00003) %>%
  na.omit()

#Convert to mm/day
df$Q<-df$Q/(gages$Watershed_Area[gages$SOURCE_FEA==gage])*(0.3048^3)/10*86400


if(nrow(df)>0){
  #Estimate baseflow using EcoHydRology package
  df$Q_bf <- BaseflowSeparation(df$Q, filter_parameter = 0.95)$bt
  
  #Estimate dQ/dt
  df<- df %>%
    mutate(dQ = Q_bf - lag(Q_bf)) %>%
    filter(dQ<0) %>%
    filter(Q_bf>0)
  
  # Take monthly means
  #df %>% group_by(month(Date)) %>% summarise(dQ=mean(dQ), Q=mean(Q)) %>% plot()  
  
  #Create model to esitmate a and b
  if(nrow(df)>0){
    model<-lm(log10(dQ*-1)~log10(Q), data=df)
    b<-model$coefficients[1]
    a<-model$coefficients[2]
    p_value<-summary(model)$coefficients[8]
  }else{
    a<- -9999
    b<- -9999
    p_value<- -9999
  }
}else{
  a<- -9999
  b<- -9999
  p_value<- -9999
}
plot(df$Q, df$dQ*-1, log="xy", 
     ylim=c(0.000001, 0.1), xlim=c(0.01, 10),
     ylab="dq/dt [mm/day]", xlab="Streamflow (mm)",
     cex.lab=14/12, cex.axis=10/12, 
     pch=19, col="#D3D3D366", main="Tuckahoe")
# points(df$Q[668:703],df$dQ[668:703]*-1, type="l", col="dark red", lwd=2)
# points(df$Q[1055:1090], df$dQ[1055:1090]*-1, type="l", lwd=2, col="dark blue")
# points(df$Q[2395:2400], df$dQ[2395:2400]*-1, type="l", lwd=2, col="orange")





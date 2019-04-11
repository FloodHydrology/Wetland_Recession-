#######################################################################
#Title: Time Series Analysis
#Coder: Nate Jones
#Date: 4/5/2019
#Description: Estimate recession coefficient for each gage
#######################################################################

#######################################################################
#Setup workspace-------------------------------------------------------
#######################################################################
#Clear workspace
remove(list=ls())

#load appropriate packages
library(dataRetrieval)
library(EcoHydRology)
library(tidyverse)

#Download data
gages<-read_csv("output/wetland_area.csv")[,2:5]
q_gage<-as_tibble(readNWISdv(siteNumbers=gages$SOURCE_FEA, parameterCd = "00060"))

#######################################################################
#Recession Analysis----------------------------------------------------
#######################################################################
#Create function to conduct recession analysis
recession_fun<-function(gage){
  #Filter to gage of interest (Choptank for now)
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
  
  #Export 
  c(gage, a, b, p_value)
}

#Execute function
output<-lapply(gages$SOURCE_FEA, recession_fun)
output<-data.frame(do.call(rbind, output))
colnames(output)<-c("SOURCE_FEA", "a", "b", "p_value")
output$a<-as.numeric(paste(output$a))
output$b<-as.numeric(paste(output$b))
output$p_value<-as.numeric(paste(output$p_value))

#put together with gages
gages<-left_join(gages, output)
gages<-gages[gages$a>0,]
gages<-gages[gages$p_value<0.001,]
gages<-na.omit(gages)

#######################################################################
#Plot------------------------------------------------------------------
#######################################################################
#Plot recession coefficients
plot((gages$Wetland_Area/gages$Watershed_Area), gages$b, 
     #Limits
     #ylim=c(0.02, 3),
     #Labels
     ylab="Recession Coefficient", xlab= "Wetland Area [%]",
     #points
     pch=19, cex.lab=14/12, cex.axis=10/12
     )

#Add regression model
model<-lm(gages$b~(gages$Wetland_Area/gages$Watershed_Area*100))
abline(model, lty=2, lwd=2, col="grey30")

#Add points
points((gages$Wetland_Area/gages$Watershed_Area*100), gages$a, pch=19)

#Add text
text(3, 1, "r=0.37\np=0.005", cex=14/12)


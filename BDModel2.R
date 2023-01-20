#==========================================================================
#Author:       Anthony Basooma
#Institution:  National Fisheries Resources Research 
#Assignment:   Biomass Dynamic Model
#====================================================================
#set working directory
setwd('E:\\NaFIRRI\\Projects\\ECOFISH\\stock assessment\\FOR LVFO')
library(deSolve)
library(extrafont)
library(dplyr)
library(ggplot2)
options(scipen = 999)

#====================================================================
#load data for cpueobs and catch
#====================================================================
catcheff <- read.csv('bdsdata.csv', fileEncoding="UTF-8-BOM") %>% 
  filter(syear<=2015) %>% mutate(cpueobs = catch/effort)

syears   <- pull(catcheff, syear)
catch    <- pull(catcheff, catch)
cpueobs  <- pull(catcheff, cpueobs)
effort  <-  pull(catcheff, effort)

#===========
#parameters
#===========
k        <-  2149478    #carrying capacity
rmax     <-  0.406976335359261       #instantaneous growth rate
q        <-  0.000000000005560441      #fishing mortality rate
p        <-  0.5
Bo       <-   k*p         # initial biomass 
params    <-  c(k, Bo, rmax, q) #param into the model

#=============

NPModel<-function(params){
  k              <-       params[1]
  Bo             <-       params[2]
  rmax           <-       params[3]
  q              <-       params[4]
  biomass        <-       Bo
  biomassvalues  <-       NULL
  cpuepredict    <-       NULL
  syears         <-       1:length(catch)
  for (syear in syears){
    bioest<-          rmax*biomass*(1-biomass/k)
    biomassvalues<-   c(biomassvalues, biomass)
    cpuepredict  <-   c(cpuepredict, q*biomass)
    biomass      <-   biomass+ bioest -catch[syear]
    biomass      <-   ifelse(biomass<0, 0, biomass)
  }
  SSE<-     sum((cpueobs-cpuepredict)^2)
  return(SSE)
}

#================
#Non linear minimization
#================
estimates <-nlm(NPModel, params,typsize = params,  iterlim=11000,
                hessian = TRUE, steptol = 1e-02) 

optim(params, NPModel)

optimize(params, NPModel)


#=========
#compute biomass with minimize biomass
k2              <-       estimates$est[1]
Bo2             <-       estimates$est[2]
rmax2           <-       estimates$est[3]
q2              <-       estimates$est[4]

biomass        <-       Bo
biomassvalues  <-       NULL
cpuepredict    <-       NULL
syears          <-      1:length(catch)

for (syear in syears){
  bioest<-     rmax*biomass*(1-biomass/k)
  biomassvalues<-   c(biomassvalues, biomass)
  cpuepredict  <-   c(cpuepredict, q*biomass)
  biomass      <-   biomass+ bioest -catch[syear]
  biomass      <-   ifelse(biomass<0, 0, biomass)
}

#============
npdata <- data.frame(syears, effort,catch, biomassvalues, cpuepredict, 
                     cpueobs) %>% 
  mutate(catchpredict = effort*cpuepredict, 
         srcpue = ((log(cpueobs)-log(cpuepredict))^2))

ggplot(npdata, aes(x= syears, y=biomassvalues))+
  geom_line()

ggplot(npdata, aes(x=syears, y=catch))+
  geom_line()

ggplot(npdata, aes(x=biomassvalues, y=cpuepredict))+
  geom_line()

#=============
##cpue
#============
ggplot(npdata, aes(x=syears, y=cpueobs))+
  geom_point()+
  geom_line(aes(x=syears, y=cpuepredict))


cor(biomassvalues, catch)

#====
#At equilbrium
biolevels <- seq(0,k,10) # biomass levels
eylevels <- rmax*biolevels*(1-biolevels/k) # equilibrium yield levels

eylevels <- ifelse(eylevels<0, 0, eylevels)

eqdata = data.frame(biolevels, eylevels)

ggplot(eqdata, aes(x=biolevels, y=eylevels))+
  geom_line()

Elevels <- seq(0,k,length.out= length(eylevels))
eldata = data.frame(Elevels, eylevels)

ggplot(eqdata, aes(x=Elevels, y=eylevels))+
  geom_line()

#==============
#1. Maximum sustainable yield 
MSY <- (rmax*k)/4 
MSY

#==================
Bmsy <- k/2
Bmsy

#===================
Emsy <- rmax/(2*q) 
Emsy

#===============
Fmsy <- rmax/2
Fmsy

#=================
Blim = k*0.25
Blim

#*****************************************************************                                                           *
#* TITLE:   Falling Creek Reservoir GLM-AED HOx submerged inflow 
#*         driver file preparation                               *
#* AUTHORS:  C.C. Carey                                          *
#* DATE:   Originally developed 16 July 2018; Last modified 11 Avril 2022                           
#* NOTES:  CCC subsequently edited on 1 June 2020 and made tidy,
#*         with subsequent tweaks to annotation in summer 2021. 
#*         SSS = side stream supersaturation, another way of describing
#*         the hypolimnetic oxygenation system (HOx)
#*****************************************************************

setwd("FCR_2013_2019GLMHistoricalRun_GLMv3beta/inputs")

#load packages
library(zoo)
library(tidyverse)
library(lubridate)
library(magrittr)

#pull in SSS operations file 
inflowoxy <- read.csv("Calc_HOX_flow_DO_20211102.csv", header=T) %>%
#within this file, the eductor nozzles increased flow rate by a factor of 4, so 227 LPM = 1135 LPM
  select(time,SSS_m3.day,mmol.O2.m3.day) %>%
  mutate(SSS_m3.day = SSS_m3.day * (1/86400))  %>%
  #convert m3/day to m3/second as needed by GLM
  mutate(SALT = rep(0,length(inflowoxy$time))) %>%
  #add salt column as needed by GLM
  mutate(time = as.POSIXct(strptime(time, "%Y-%m-%d", tz="EST"))) %>%
  rename(FLOW = SSS_m3.day, OXY_oxy=mmol.O2.m3.day)

#some diagnostic plots
plot(inflowoxy$time, inflowoxy$FLOW, type = "o")
plot(inflowoxy$time, inflowoxy$OXY_oxy, type = "l", col = "red", ylab="mmol O2/m3/d added per day",
     xlab="Time")

#now, need to inject high levels of O2 that match HOx operations, with low flow rate and no solutes
highox <- inflowoxy 
highox <- inflowoxy %>% 
  mutate(load = FLOW*OXY_oxy) %>% 
  mutate(FLOW_lo = rep(1e-08, length(highox$FLOW))) %>% #using flow rate of 1e-08
  mutate(OXY_oxy_lo = load/FLOW_lo) %>% 
  select(time, FLOW_lo, SALT, OXY_oxy_lo) %>% 
  rename(FLOW = FLOW_lo, OXY_oxy = OXY_oxy_lo) %>% 
  mutate(NIT_amm = rep(0, length(highox$time)),
         NIT_nit = rep(0, length(highox$time)),
         PHS_frp = rep(0, length(highox$time)),
         OGM_doc = rep(0, length(highox$time)),
         OGM_docr = rep(0, length(highox$time)),
         OGM_poc = rep(0, length(highox$time)),
         OGM_don = rep(0, length(highox$time)),
         OGM_donr = rep(0, length(highox$time)),
         OGM_pon = rep(0, length(highox$time)),
         OGM_dop = rep(0, length(highox$time)),
         OGM_dopr = rep(0, length(highox$time)),
         OGM_pop = rep(0, length(highox$time)),
         CAR_dic = rep(0, length(highox$time)),
         CAR_pH = rep(0, length(highox$time)),
         CAR_ch4 = rep(0, length(highox$time)),
         SIL_rsi = rep(0, length(highox$time)),
         TRC_tr1 = rep(1,length(highox$time)),
         TRC_age = rep(0,length(highox$time)),
         NCS_ss1 = rep(0,length(highox$time)),
         NCS_ss2 = rep(0,length(highox$time)),
         CAR_ch4_bub = rep(0,length(highox$time)),
         PHS_frp_ads = rep(0,length(highox$time)),
         OGM_cpom = rep(0,length(highox$time)),
         PHY_cyano = rep(0,length(highox$time)),
         PHY_cyano_IN = rep(0,length(highox$time)),
         PHY_cyano_IP = rep(0,length(highox$time)),
         PHY_green = rep(0,length(highox$time)),
         PHY_green_IN = rep(0,length(highox$time)),
         PHY_green_IP = rep(0,length(highox$time)),
         PHY_diatom = rep(0,length(highox$time)),
         PHY_diatom_IN = rep(0,length(highox$time)),
         PHY_diatom_IP = rep(0,length(highox$time)),
         BIV_filtfrac = rep(0,length(highox$time)))

#now, need to set water temperature of this file to CTD observations at 8 m, the depth
# the HOx injects water into the hypolimnion
#inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/12/0a62d1946e8d9a511bc1404e69e59b8c" 
#infile1 <- paste0(getwd(),"/CTD_2013_2021.csv")
#download.file(inUrl1,infile1,method="curl")

CTD<-read.csv("CTD_2013_2021.csv", header=TRUE) #now need to get temp at 8m for inflow
CTD8 <- CTD %>%
  select(Reservoir:Temp_C) %>%
  dplyr::filter(Reservoir=="FCR") %>%
  dplyr::filter(Site==50) %>%
  rename(time=Date, TEMP=Temp_C) %>%
  mutate(time = as.POSIXct(strptime(time, "%Y-%m-%d", tz="EST"))) %>%
  mutate(Depth_m = round(Depth_m, digits=0)) %>%
  group_by(time) %>% 
  dplyr::filter(Depth_m==8) %>%
  summarise(TEMP=mean(TEMP)) %>%
  dplyr::filter(TEMP<17) #remove outlier from 2014

#diagnostic plot to check for 8m water temp
plot(CTD8$time, CTD8$TEMP, type = "o")

#make final SSS inflow file, setting Dec 31 of each year to 4oC in lieu of CTD data for interpolation
SSS_inflowALL<-merge(highox,CTD8, by="time",all.x=TRUE)
SSS_inflowALL$TEMP[231]<-4
SSS_inflowALL$TEMP[596]<-4
SSS_inflowALL$TEMP[961]<-4
SSS_inflowALL$TEMP[1326]<-4
SSS_inflowALL$TEMP[1691]<-4
SSS_inflowALL$TEMP[2056]<-4
SSS_inflowALL$TEMP[2422]<-4 
SSS_inflowALL$TEMP[2788]<-4 
SSS_inflowALL$TEMP[3153]<-4 #set last row as 4oC in prep for freezing
SSS_inflowALL$TEMP<-na.fill(na.approx(SSS_inflowALL$TEMP),"extend")
plot(SSS_inflowALL$time, SSS_inflowALL$TEMP, type = "o")

#get everything in order
SSS_inflowALL1<-SSS_inflowALL %>%
  select(time, FLOW, TEMP, SALT, 
         #TRC_tr1, TRC_age, 
         #NCS_ss1, 
         #NCS_ss2, 
         OXY_oxy, CAR_dic, 
         CAR_pH, CAR_ch4, #CAR_ch4_bub, 
         SIL_rsi, NIT_amm, NIT_nit, PHS_frp, 
         #PHS_frp_ads,
         OGM_doc, OGM_poc, OGM_don, OGM_pon, OGM_dop, OGM_pop, OGM_docr, OGM_donr, 
         OGM_dopr, OGM_cpom, 
         PHY_cyano, PHY_cyano_IN, PHY_cyano_IP, PHY_green, 
         PHY_green_IN, PHY_green_IP, PHY_diatom, PHY_diatom_IN, PHY_diatom_IP
         #, BIV_filtfrac
         ) 
#get all of the columns in order

SSS_inflowALL1[which(duplicated(SSS_inflowALL1$time)),] #identify if there are repeated dates
SSS_inflowALL1 <- SSS_inflowALL1[(!duplicated(SSS_inflowALL1$time)),] #remove repeated dates

#et voila! the final observed inflow file for the SSS for 2 pools of DOC
write.csv(SSS_inflowALL1, "FCR_SSS_inflow_2013_2021_20220413_allfractions_2DOCpools.csv", row.names = FALSE)

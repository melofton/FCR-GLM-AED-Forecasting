#*****************************************************************                                                           *
#* TITLE:   Falling Creek Reservoir GLM-AED climatology stream  
#*          inflow file preparation                                          
#* AUTHORS:  C.C. Carey                                          
#* DATE:   Originally developed 10 Nov 2021; Last modified 12 Nov 2021                            
#* NOTES:  CCC developed to estimate reservoir inflows for FCR; 
#*         then converted over for making the climatology inflow file
#*****************************************************************

setwd("FCR_2013_2019GLMHistoricalRun_GLMv3beta/inputs")
setwd("./inputs")

if(!require('pacman')) install.packages('pacman'); library('pacman')
pacman::p_load(zoo, EcoHydRology, rMR, tidyverse, magrittr, lubridate, 
               ncdf4, glmtools, GLM3r, birk, dplyr)
library(tidyverse)

#set the start and end dates of the data going into potential aggregation for
# your inflow file (you'll nail down exact dates below)
start_date<-as.POSIXct(strptime("2015-07-07", "%Y-%m-%d", tz="EST"))
end_date<-as.POSIXct(strptime("2020-12-31", "%Y-%m-%d", tz="EST"))

#creating new dataframe with list of all dates
datelist<-seq.Date(as.Date(start_date),as.Date(end_date), "days") #changed from May 15, 2013 because of NA in flow
datelist<-as.data.frame(datelist)
colnames(datelist)=c("time")
datelist$time<-as.POSIXct(strptime(datelist$time, "%Y-%m-%d", tz="EST"))

#first read in FCR weir inflow file from EDI (updated for 2013-Dec 2020)
if(!file.exists('inflow_for_EDI_2013_10Jan2021.csv')){
  inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/202/7/f5fa5de4b49bae8373f6e7c1773b026e" 
  infile1 <- paste0(getwd(),"/inflow_for_EDI_2013_10Jan2021.csv")
  download.file(inUrl1,infile1,method="curl")  
}

inflow<-readr::read_csv("inflow_for_EDI_2013_10Jan2021.csv") %>% 
  dplyr::select(DateTime, WVWA_Flow_cms, WVWA_Temp_C) %>% 
  dplyr::rename(time=DateTime, FLOW=WVWA_Flow_cms, TEMP=WVWA_Temp_C) %>%
  mutate(time = as.POSIXct(strptime(time, "%Y-%m-%d", tz="EST"))) %>%
  dplyr::filter(time < end_date + lubridate::days(1)) %>%
  group_by(time) %>% 
  dplyr::summarise(FLOW=mean(FLOW), TEMP=mean(TEMP)) #gives averaged daily flow per day in m3/s

#diagnostic plot
# plot(inflow$time, inflow$FLOW)

#merge inflow file with datelist to make sure that we have all days covered 
#interpolating the few missing days
weir <- merge(datelist, inflow, by="time", all.x=TRUE) %>%
  mutate(FLOW = na.fill(na.approx(FLOW),"extend")) %>%
  mutate(TEMP = na.fill(na.approx(TEMP),"extend")) %>%
  mutate(SALT = rep(0,length(datelist)))

#some diagnostic plots of inflow weir
# plot(weir$time, weir$FLOW, type = "o")
#par(new = T)
# plot(weir$time, weir$TEMP, type = "l", col = "red")

#now let's merge with chemistry
#first pull in FCR chem data from 2013-2020 from EDI
#pull in FCR chem data from 2013-2020 from EDI
if(!file.exists("chemistry_2013_2020.csv")){
  inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/9/fe500aac19d1a0d78bb2cb1d196cdbd7" 
  infile1 <- paste0(getwd(),"/chemistry_2013_2020.csv")
  download.file(inUrl1,infile1,method="curl")
}

FCRchem <- readr::read_csv("chemistry_2013_2020.csv") %>% 
  select(Reservoir:DIC_mgL) %>%
  dplyr::filter(Reservoir=="FCR") %>%
  dplyr::filter(Site==100) %>% #inflow site code
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>%
  rename(time = DateTime) %>%
  dplyr::filter(TP_ugL < 100) %>% #remove outliers
  select(time:DIC_mgL) 

#read in lab dataset of dissolved silica, measured in summer 2014 only
if(!file.exists("silica_master_df.csv")){
  inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/542/1/791ec9ca0f1cb9361fa6a03fae8dfc95" 
  infile1 <- paste0(getwd(),"/silica_master_df.csv")
  download.file(inUrl1,infile1,method="curl")
}

silica <- read.csv("silica_master_df.csv", header=T) %>%
  dplyr::filter(Reservoir == "FCR") %>% 
  dplyr::filter(Site == 100) %>% #100 = weir inflow site
  select(DateTime, DRSI_mgL) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>%
  rename(time = DateTime)

#diagnostic plot of silica
# plot(silica$time, silica$DRSI_mgL)
# hist(silica$DRSI_mgL)
# median(silica$DRSI_mgL) 
#median conc is going to be set as the constant Si inflow conc in inflow

all_data<-merge(weir, FCRchem, by="time", all.x=TRUE)

#read in dataset of CH4 from EDI
if(!file.exists("Dissolved_CO2_CH4_Virginia_Reservoirs.csv")){
  inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/551/5/38d72673295864956cccd6bbba99a1a3" 
  infile1 <- paste0(getwd(),"/Dissolved_CO2_CH4_Virginia_Reservoirs.csv")
  download.file(inUrl1,infile1,method="curl")
}

ghg <- read.csv("Dissolved_CO2_CH4_Virginia_Reservoirs.csv", header=T) %>%
  dplyr::filter(Reservoir == "FCR") %>%
  dplyr::filter(Site == 100) %>% #weir inflow
  select(DateTime, ch4_umolL) %>%
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>%
  rename(time = DateTime, CAR_ch4 = ch4_umolL) %>%
  group_by(time) %>%
  drop_na %>% 
  summarise(CAR_ch4 = mean(CAR_ch4)) %>%
  dplyr::filter(CAR_ch4<0.2) #remove outliers
# plot(ghg$time, ghg$CAR_ch4)

ghg1 <- merge(datelist, ghg, by="time", all.x=TRUE) 
#need to interpolate missing data, but first need to fill first & last values
ghg1$CAR_ch4[1]<-ghg$CAR_ch4[which.closest(ghg$time, start_date)]
ghg1$CAR_ch4[length(ghg1$CAR_ch4)]<-ghg$CAR_ch4[which.closest(ghg$time, end_date)]
ghg1$CAR_ch4 <- na.fill(na.approx(ghg1$CAR_ch4), "extend")
# plot(ghg1$time, ghg1$CAR_ch4) 

#some other cool long-term plots
# plot(alldata$time, alldata$SRP_ugL)
# plot(alldata$time, alldata$DOC_mgL)
# plot(alldata$time, alldata$NO3NO2_ugL) #something's off here; need Heather's modified EDI chem file
# plot(alldata$time, alldata$NH4_ugL)
# plot(alldata$time, alldata$TN_ugL)
# plot(alldata$time, alldata$TP_ugL)
# plot(alldata$time, alldata$DIC_mgL)

alldata<-merge(all_data, ghg1, by="time", all.y=TRUE) %>% 
  group_by(time) %>% 
  summarise_all(mean, na.RM=TRUE)
#merge chem with CH4 data, truncating to start and end date period of CH4 (not chem)

lastrow <- length(alldata$time) #need for extend function below
#now need to interpolate missing values in chem; setting 1st and last value in time series as medians
#then linearly interpolating the middle missing values
alldata$TN_ugL[1]<-median(na.exclude(alldata$TN_ugL))
alldata$TN_ugL[lastrow]<-median(na.exclude(alldata$TN_ugL)) 
alldata$TN_ugL<-na.fill(na.approx(alldata$TN_ugL),"extend")

alldata$TP_ugL[1]<-median(na.exclude(alldata$TP_ugL))
alldata$TP_ugL[lastrow]<-median(na.exclude(alldata$TP_ugL))
alldata$TP_ugL<-na.fill(na.approx(alldata$TP_ugL),"extend")

alldata$NH4_ugL[1]<-median(na.exclude(alldata$NH4_ugL))
alldata$NH4_ugL[lastrow]<-median(na.exclude(alldata$NH4_ugL))
alldata$NH4_ugL<-na.fill(na.approx(alldata$NH4_ugL),"extend")

alldata$NO3NO2_ugL[1]<-median(na.exclude(alldata$NO3NO2_ugL))
alldata$NO3NO2_ugL[lastrow]<-median(na.exclude(alldata$NO3NO2_ugL))
alldata$NO3NO2_ugL<-na.fill(na.approx(alldata$NO3NO2_ugL),"extend")

alldata$SRP_ugL[1]<-median(na.exclude(alldata$SRP_ugL))
alldata$SRP_ugL[lastrow]<-median(na.exclude(alldata$SRP_ugL))
alldata$SRP_ugL<-na.fill(na.approx(alldata$SRP_ugL),"extend")

alldata$DOC_mgL[1]<-median(na.exclude(alldata$DOC_mgL))
alldata$DOC_mgL[lastrow]<-median(na.exclude(alldata$DOC_mgL))
alldata$DOC_mgL<-na.fill(na.approx(alldata$DOC_mgL),"extend")

alldata$DIC_mgL[1]<-median(na.exclude(alldata$DIC_mgL))
alldata$DIC_mgL[lastrow]<-median(na.exclude(alldata$DIC_mgL))
alldata$DIC_mgL<-na.fill(na.approx(alldata$DIC_mgL),"extend")

alldata <- alldata[(!duplicated(alldata$time)),]#remove duplicated dates

#need to convert mass observed data into mmol/m3 units for two pools of organic carbon
weir_inflow <- alldata %>% 
  mutate(NIT_amm = NH4_ugL*1000*0.001*(1/18.04)) %>% 
  mutate(NIT_nit = NO3NO2_ugL*1000*0.001*(1/62.00)) %>% #as all NO2 is converted to NO3
  mutate(PHS_frp = SRP_ugL*1000*0.001*(1/94.9714)) %>% 
  mutate(OGM_doc = DOC_mgL*1000*(1/12.01)* 0.10) %>% #assuming 10% of total DOC is in labile DOC pool (Wetzel page 753)
  mutate(OGM_docr = 1.5*DOC_mgL*1000*(1/12.01)* 0.90) %>% #assuming 90% of total DOC is in recalcitrant DOC pool
  mutate(TN_ugL = TN_ugL*1000*0.001*(1/14)) %>% 
  mutate(TP_ugL = TP_ugL*1000*0.001*(1/30.97)) %>% 
  mutate(OGM_poc = 0.1*(OGM_doc+OGM_docr)) %>% #assuming that 10% of DOC is POC (Wetzel page 755)
  mutate(OGM_don = (5/6)*(TN_ugL-(NIT_amm+NIT_nit))*0.10) %>% #DON is ~5x greater than PON (Wetzel page 220)
  mutate(OGM_donr = (5/6)*(TN_ugL-(NIT_amm+NIT_nit))*0.90) %>% #to keep mass balance with DOC, DONr is 90% of total DON
  mutate(OGM_pon = (1/6)*(TN_ugL-(NIT_amm+NIT_nit))) %>% #detemined by subtraction
  mutate(OGM_dop = 0.3*(TP_ugL-PHS_frp)*0.10) %>% #Wetzel page 241, 70% of total organic P = particulate organic; 30% = dissolved organic P
  mutate(OGM_dopr = 0.3*(TP_ugL-PHS_frp)*0.90) %>% #to keep mass balance with DOC & DON, DOPr is 90% of total DOP
  mutate(OGM_pop = 10*TP_ugL) %>% # #In lieu of having the adsorbed P pool activated in the model, need to have higher complexed P
  mutate(CAR_dic = DIC_mgL*1000*(1/52.515)) #Long-term median pH of FCR is 6.5, at which point CO2/HCO3 is about 50-50
#given this disparity, using a 50-50 weighted molecular weight (44.01 g/mol and 61.02 g/mol, respectively)


#reality check of mass balance
# hist(weir_inflow$TP_ugL - (weir_inflow$PHS_frp + weir_inflow$OGM_dop + weir_inflow$OGM_dopr + weir_inflow$OGM_pop))
# hist(weir_inflow$TN_ugL - (weir_inflow$NIT_amm + weir_inflow$NIT_nit + weir_inflow$OGM_don + weir_inflow$OGM_donr + weir_inflow$OGM_pon))

#creating OXY_oxy column using RMR package, assuming that oxygen is at 100% saturation in this very well-mixed stream
for(i in 1:length(weir_inflow$TEMP)){
  weir_inflow$OXY_oxy[i]<-(temp.C= Eq.Ox.conc(weir_inflow$TEMP[i], elevation.m = 506,
                                              bar.press = NULL, bar.units = NULL,
                                              out.DO.meas = "mg/L",
                                              salinity = 0, salinity.units = "pp.thou"))*1000*(1/32)
}

#clean it up and get vars in order
weir_inflow <- weir_inflow %>%
  select(time, FLOW, TEMP, SALT, OXY_oxy, NIT_amm:CAR_dic, CAR_ch4) %>% 
  mutate(SIL_rsi = rep(median(silica$DRSI_mgL),length(weir_inflow$time))) %>%
  mutate(SIL_rsi = SIL_rsi*1000*(1/60.08)) %>% #setting the Silica concentration to the median 2014 inflow concentration for consistency
  mutate_if(is.numeric, round, 4) #round to 4 digits 


##############################################
#First, let's make an outflow for weir-only outflow
outflow <- weir_inflow %>% #from above: this has both stream inflows together
  select(time, FLOW) %>%
  mutate_if(is.numeric, round, 4) #round to 4 digits 

#diagnostic plot
# plot(outflow$time, outflow$FLOW)

#write file
write.csv(outflow, "FCR_spillway_outflow_WeirOnly_2015_2020_20211114.csv", row.names=F)
##############################################

#now, let's create the climatology weir inflow file for Quinn, by taking the
# average of day of year for all analytes
# First, need to choose which data go into the climatology calculations

clima_start <- as.POSIXct(strptime("2015-07-07", "%Y-%m-%d", tz="EST"))
clima_end <-as.POSIXct(strptime("2018-07-08", "%Y-%m-%d", tz="EST"))

weir_inflow_dates <- weir_inflow %>% 
  dplyr::filter(time>clima_start & time<clima_end) %>% 
  mutate(DOY = yday(time)) %>%
  select(time, DOY) 

meas_vars <- weir_inflow %>% 
  select(time:OXY_oxy)

#check that all vars looking ok
# plot(weir_inflow$time, weir_inflow$FLOW)
# plot(weir_inflow$time, weir_inflow$TEMP)
# plot(weir_inflow$time, weir_inflow$SALT)
# plot(weir_inflow$time, weir_inflow$OXY_oxy)
# plot(weir_inflow$time, weir_inflow$NIT_amm)#flagged only after 2015
# plot(weir_inflow$time, weir_inflow$NIT_nit)#flagged only after 2015
# plot(weir_inflow$time, weir_inflow$PHS_frp)#flagged only after 2015
# plot(weir_inflow$time, weir_inflow$OGM_doc)#flagged only after 2015
# plot(weir_inflow$time, weir_inflow$OGM_docr)#flagged only after 2015
# plot(weir_inflow$time, weir_inflow$OGM_poc)#flagged only after 2015
# plot(weir_inflow$time, weir_inflow$OGM_don)#flagged only after 2015
# plot(weir_inflow$time, weir_inflow$OGM_donr)#flagged only after 2015
# plot(weir_inflow$time, weir_inflow$OGM_pon)#flagged only after 2015
# plot(weir_inflow$time, weir_inflow$OGM_dop)##flagged only after 2015
# plot(weir_inflow$time, weir_inflow$OGM_pop)##flagged only after 2015
# plot(weir_inflow$time, weir_inflow$CAR_dic)##flagged only after 2015
# plot(weir_inflow$time, weir_inflow$CAR_ch4)#flagged only after 2015
# plot(weir_inflow$time, weir_inflow$SIL_rsi)

#now make mean climatology
mean_DOY_data <- weir_inflow %>% 
  dplyr::filter(time> clima_start & time< clima_end) %>% 
  mutate(DOY = yday(time)) %>%
  group_by(DOY) %>%
  dplyr::summarise(across(c("NIT_amm":"SIL_rsi"),mean))

climatology_mean <- dplyr::left_join(weir_inflow_dates, mean_DOY_data) 
climatology_mean1 <- merge(meas_vars, climatology_mean, by="time") %>% 
  select(!(DOY))

climatology_mean2 <- climatology_mean1 %>% 
  mutate(OGM_docr1 = -238.5586 + 4101.3976*FLOW + 2.1472*OGM_docr + (-19.1272*OGM_docr*FLOW)) %>% 
  select(time:OGM_doc,OGM_docr1,OGM_poc:SIL_rsi) %>% 
  rename(OGM_docr=OGM_docr1)
#this is my model, where I predict what stream OGM docr concentrations need to
# be based off of my 

write.csv(climatology_mean2, "climatology_mean_weir_inflow_higherDOCr.csv",row.names = F)







######## code that isn't needed now!
#climatology median inflow
median_DOY_data <- weir_inflow %>% 
  dplyr::filter(time > clima_start & time < clima_end) %>% 
  mutate(DOY = yday(time)) %>%
  group_by(DOY) %>%
  dplyr::summarise(across(c("NIT_amm":"SIL_rsi"),median))

climatology_median <- dplyr::left_join(weir_inflow_dates, median_DOY_data) 
climatology_median1 <- merge(meas_vars, climatology_median, by="time") %>% 
  select(!(DOY))

climatology_mean2 <- climatology_mean1 %>% 
  mutate(OGM_docr1 = -238.5586 + 4101.3976*FLOW + 2.1472*OGM_docr + (-19.1272*OGM_docr*FLOW)) %>% 
  select(FLOW:OGM_doc,OGM_docr1,OGM_poc:SIL_rsi) %>% 
  rename(OGM_docr=OGM_docr1)

write.csv(climatology_median1, "climatology_median_weir_inflow.csv",row.names = F)
######## 







# #############################################
# extra code for modifying stream inflows
# #modifying inflow DOCr concentration to better predict 1.6m lake concentrations
# 
# weir_inflow1 <- weir_inflow %>% 
#   mutate(OGM_docr1 = -238.5586 + 4101.3976*FLOW + 2.1472*OGM_docr + (-19.1272*OGM_docr*FLOW))
# 
# Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)               -238.5586    57.4591  -4.152 7.14e-05 ***
#   comp2$FLOW                4101.3976  1577.4918   2.600  0.01080 *  
#   comp2$OGM_docr               2.1472     0.2512   8.548 1.95e-13 ***
#   comp2$FLOW:comp2$OGM_docr  -19.1272     6.2691  -3.051  0.00295 ** 
#   
# 
# ######building a model to compare stream inflow DOC concentrations with lake data
# #inlake 1.6 concentrations
# obs<-read.csv('./field_data/field_chem.csv', header=TRUE) %>% #read in observed chemistry data
#   dplyr::mutate(time = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>%
#   dplyr::filter(Depth==1.6) %>% 
#   select(time, OGM_docr)
# 
# #stream concs
# weirgrab<-weir_inflow %>% 
#   select(time,FLOW, TEMP, OGM_docr)
# 
# #modeled
# nc_file <- file.path("/Users/cayelan/Dropbox/ComputerFiles/SCC/AED_Forecasting/FCR-GLM-AED-Forecasting/FCR_2013_2019GLMHistoricalRun_GLMv3beta/output", 'output.nc') #defines the output.nc file 
# 
# mod<- get_var(nc_file, "OGM_docr", reference="surface", z_out=1.6) %>%
#   mutate(time = as.POSIXct(strptime(DateTime, "%Y-%m-%d", tz="EST"))) %>% 
#   select(time,OGM_docr_1.6)
# 
# comp<-base::merge(obs,weirgrab,by=c("time"))
# comp1<-base::merge(comp,mod, by=c("time")) %>% 
#   rename(lakeobs=OGM_docr.x, 
#          stream=OGM_docr.y,
#          mod=OGM_docr_1.6) %>% 
#   mutate(factor=lakeobs-mod)
# 
# climatology <- climatology_mean1 %>% 
#   select(time, OGM_docr) %>% 
#   mutate(DOY = yday(time))
# 
# plot(climatology$time, climatology$OGM_docr)
# 
# comp2<-base::merge(comp1,climatology, by=c("time"))
# 
# plot(comp2$time, comp2$lakeobs, col="red", type="o")
# lines(comp2$time, comp2$stream, col="blue")
# lines(comp2$time, comp2$mod, col="black")
# plot(comp2$time,comp2$factor)
# plot(comp2$time, comp2$OGM_docr)
# 
# model<- lm(comp2$stream ~ comp2$OGM_docr*comp2$FLOW)
# summary(model)
# 
# model<- lm(comp2$stream ~ comp2$FLOW*comp2$OGM_docr)
# 
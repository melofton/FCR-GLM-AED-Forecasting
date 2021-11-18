#*****************************************************************                                                           *
#* TITLE:   Falling Creek Reservoir GLM-AED forecasting water quality  
#*          variable conversions                                    
#* AUTHORS: C.C. Carey                                          
#* DATE:   Originally developed 16 Nov 2021                            
#* NOTES:  CCC developed to estimate compare observed DOC with fDOM EXO,
#*         chla CTD --> EXO, DO CTD --> DO EXO
#*****************************************************************

setwd("FCR_2013_2019GLMHistoricalRun_GLMv3beta/inputs") #OR
setwd("./inputs")

if(!require('pacman')) install.packages('pacman'); library('pacman')
pacman::p_load(zoo, EcoHydRology, rMR, tidyverse, magrittr, lubridate, 
               ncdf4, glmtools, GLM3r, readr, dplyr)

#first, we need to convert fDOM into DOC total units
#####################################################
#pull in FCR chem data from 2013-2020 from EDI
if(!file.exists("chemistry_2013_2020.csv")){
  inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/9/fe500aac19d1a0d78bb2cb1d196cdbd7" 
  infile1 <- paste0(getwd(),"/chemistry_2013_2020.csv")
  download.file(inUrl1,infile1,method="curl")
}

FCRchem <- readr::read_csv("chemistry_2013_2020.csv") %>%
  dplyr::filter(Reservoir == "FCR", 
                Site == 50,  #deep hole site code
                Depth_m == 1.6,
                DOC_mgL < 10) %>% 
  dplyr::rename(time = DateTime)  %>% 
  dplyr::mutate(time = force_tz(time, "EST"),
                time = with_tz(time, "UTC"),
                time = date(time)) %>%
  dplyr::select(time, DOC_mgL) %>% 
  na.exclude()

#pull in catwalk sensor data from 2018-2020 from EDI
if(!file.exists("Catwalk_EDI_2020.csv")){
  inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/271/5/c1b1f16b8e3edbbff15444824b65fe8f" 
  infile1 <- paste0(getwd(),"/Catwalk_EDI_2020.csv")
  download.file(inUrl1,infile1,method="curl")
}

catwalk <- readr::read_csv("Catwalk_EDI_2020.csv") %>%
  dplyr::filter(Reservoir == "FCR", 
                Site == 50) %>%  #deep hole site code
  dplyr::rename(time = DateTime)  %>% 
  dplyr::mutate(time = force_tz(time, "EST"),
                time = with_tz(time, "UTC")) %>%
  dplyr::select(time, EXODO_mgL_1, RDO_mgL_5_adjusted, RDO_mgL_9_adjusted, 
                EXOChla_ugL_1, EXOfDOM_QSU_1) %>% 
  dplyr::mutate(time = date(time)) %>% 
  dplyr::group_by(time) %>% 
  dplyr::summarise(across(c("EXODO_mgL_1":"EXOfDOM_QSU_1"),mean)) %>% 
  dplyr::filter(!(time %in% as_date(c("2018-10-08","2018-10-22","2019-10-15", 	
                                      "2020-06-03","2020-10-28")))) #removing outliers & turnover data

fdom2doc <- dplyr::left_join(catwalk, FCRchem, by = c('time')) %>% 
  select(time, EXOfDOM_QSU_1, DOC_mgL) %>% 
  mutate(OGM_doc_total = DOC_mgL*1000*(1/12.01)) %>% #conversion factor to get to mmoles/m3
  na.exclude()

fdom2doc_model <- lm(OGM_doc_total ~ EXOfDOM_QSU_1, data = fdom2doc)
#model from which summary stats will be extracted

#figure of fDOM vs OGM_doc_total
ggplot(fdom2doc, aes(EXOfDOM_QSU_1, OGM_doc_total, colour=factor(time))) + 
  geom_point(alpha = 0.6) + 
  geom_abline(intercept = as.numeric(fdom2doc_model$coefficients[1]),
              slope = as.numeric(fdom2doc_model$coefficients[2]))

#equation for the fDOM EXO --> OGM_doc_total model
print(fdom2doc_model$coefficients[1])#intercept
print(fdom2doc_model$coefficients[2])#slope
print(sd(fdom2doc_model$residuals)) #standard deviation of the residuals


# #the same regression, except using DOC mg/L instead of OGM_doc_total in molar units
# fdom2doc_model <- lm(DOC_mgL ~ EXOfDOM_QSU_1, data = fdom2doc)
# 
# ggplot(fdom2doc, aes(EXOfDOM_QSU_1, DOC_mgL, colour=factor(time))) + 
#   geom_point(alpha = 0.6) + 
#   geom_abline(intercept = as.numeric(fdom2doc_model$coefficients[1]),
#               slope = as.numeric(fdom2doc_model$coefficients[2]), colour = "black") +
#   geom_abline(intercept = -0.4675,
#               slope = 0.2575, colour = "red")
#black = our regression
#red = Howard et al. 2021 AS


#####################################################
#convert 1.6 EXO chla into CTD chla

#need to pull in CTD data from EDI
if(!file.exists("CTD_final_2013_2020.csv")){
  inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/11/d771f5e9956304424c3bc0a39298a5ce" 
  infile1 <- paste0(getwd(),"/CTD_final_2013_2020.csv")
  download.file(inUrl1,infile1,method="curl")
}

#extracting DO and chla CTD for further analysis data 
ctd <- read_csv("CTD_final_2013_2020.csv") %>% 
  dplyr::filter(Reservoir == "FCR", 
                Site == 50) %>%
  dplyr::select(Date, Depth_m, DO_mgL, Chla_ugL) %>% 
  dplyr::rename(time = Date, depth = Depth_m) %>% 
  dplyr::mutate(time = force_tz(time, "EST"),
                time = with_tz(time, "UTC"))

#need to filter CTD data at EXO depth (1.6 m)
ctd_1.6 <- ctd %>% 
  dplyr::filter(depth > 1.5,
                depth < 1.7) %>% 
  dplyr::mutate(time = date(time)) %>% 
  dplyr::group_by(time) %>% 
  dplyr::summarise(across(c("DO_mgL":"Chla_ugL"),mean))

EXO2ctd <- dplyr::inner_join(ctd_1.6, catwalk, by = c('time'))

#let's do chla conversion of 1.6 EXO chla into CTD chla first
EXO2ctd_chla <- EXO2ctd %>% 
  dplyr::select(time, Chla_ugL, EXOChla_ugL_1) %>% 
  na.exclude()

#model from which summary stats will be extracted
EXO2ctd_chla_model <- lm(Chla_ugL ~ EXOChla_ugL_1, data=EXO2ctd_chla)
  
#let's look at that chla!
ggplot(EXO2ctd_chla, aes(EXOChla_ugL_1, Chla_ugL, colour=factor(time))) + 
  geom_point(alpha = 0.6) + 
  geom_abline(intercept = as.numeric(EXO2ctd_chla_model$coefficients[1]),
              slope = as.numeric(EXO2ctd_chla_model$coefficients[2]))

#summary stats for chla model
print(EXO2ctd_chla_model$coefficients[1])#intercept
print(EXO2ctd_chla_model$coefficients[2])#slope
print(sd(EXO2ctd_chla_model$residuals)) #standard deviation of the residuals

#####################################################
#convert 1.6 m EXO DO to CTD DO

#let's do DO conversion of 1.6 EXO DO into CTD DO second
EXO2ctd_do <- EXO2ctd %>% 
  dplyr::select(time, DO_mgL, EXODO_mgL_1) %>% 
  na.exclude() %>% 
  dplyr::mutate(OXY_oxy_CTD = DO_mgL*1000*(1/32),
                OXY_oxy_EXO = EXODO_mgL_1*1000*(1/32))

#model from which summary stats will be extracted
EXO2ctd_do_model <- lm(OXY_oxy_CTD ~ OXY_oxy_EXO, data=EXO2ctd_do)

#let's look at that DO!
ggplot(EXO2ctd_do, aes(OXY_oxy_EXO, OXY_oxy_CTD, colour=factor(time))) + 
  geom_point(alpha = 0.6) + 
  geom_abline(intercept = as.numeric(EXO2ctd_do_model$coefficients[1]),
              slope = as.numeric(EXO2ctd_do_model$coefficients[2]))

#summary stats for DO model
print(EXO2ctd_do_model$coefficients[1])#intercept
print(EXO2ctd_do_model$coefficients[2])#slope
print(sd(EXO2ctd_do_model$residuals)) #standard deviation of the residuals
#THIS REGRESSION IS SO BAD THAT WE'RE NOT GOING TO USE IT IN FORECASTING MODE
# BECAUSE DO IS SO UNDERPREDICTED BY CTD... GOING TO PRETEND EXO IS REAL-LIFE 
# WHICH IS PROBABLY CLOSER TO REALITY, ANYWAY!

#####################################################
#convert 5 m DO sensor to CTD DO

#need to filter CTD data at sensor depth (5 m)
ctd_5 <- ctd %>% 
  dplyr::filter(depth > 4.9,
                depth < 5.1) %>% 
  dplyr::mutate(time = date(time)) %>% 
  dplyr::group_by(time) %>% 
  dplyr::summarise(CTD_DO_5 = mean(DO_mgL))

EXO2ctd_5 <- dplyr::inner_join(ctd_5, catwalk, by = c('time'))

#let's do  conversion of 5 m do into CTD do first
EXO2ctd_do5 <- EXO2ctd_5 %>% 
  dplyr::select(time, CTD_DO_5, RDO_mgL_5_adjusted) %>% 
  na.exclude()

#model from which summary stats will be extracted
EXO2ctd_do_5_model <- lm(CTD_DO_5 ~ RDO_mgL_5_adjusted, data=EXO2ctd_do5)

#let's look at that chla!
ggplot(EXO2ctd_do5, aes(RDO_mgL_5_adjusted, CTD_DO_5, colour=factor(time))) + 
  geom_point(alpha = 0.6) + 
  geom_abline(intercept = as.numeric(EXO2ctd_do_5_model$coefficients[1]),
              slope = as.numeric(EXO2ctd_do_5_model$coefficients[2]), colour="black") +
  geom_abline(intercept = 0, slope = 1, colour = "red")#1:1 line

#summary stats for chla model
print(EXO2ctd_do_5_model$coefficients[1])#intercept
print(EXO2ctd_do_5_model$coefficients[2])#slope
print(sd(EXO2ctd_do_5_model$residuals)) #standard deviation of the residuals
#THE REGRESSION MODEL AND ONE-ONE LINE ARE SO CLOSE, NO NEED TO CONVERT EXO --> CTD AT 5M!


#####################################################
#convert 9 m DO sensor to CTD DO

#need to filter CTD data at sensor depth (9 m)
ctd_9 <- ctd %>% 
  dplyr::filter(depth > 8.9,
                depth < 9.1) %>% 
  dplyr::mutate(time = date(time)) %>% 
  dplyr::group_by(time) %>% 
  dplyr::summarise(CTD_DO_9 = mean(DO_mgL))

EXO2ctd_9 <- dplyr::inner_join(ctd_9, catwalk, by = c('time'))

#let's do do conversion of 9 m DO  into CTD do first
EXO2ctd_do9 <- EXO2ctd_9 %>% 
  dplyr::select(time, CTD_DO_9, RDO_mgL_9_adjusted) %>% 
  na.exclude()

#model from which summary stats will be extracted
EXO2ctd_do_9_model <- lm(CTD_DO_9 ~ RDO_mgL_9_adjusted, data=EXO2ctd_do9)

#let's look at that chla!
ggplot(EXO2ctd_do9, aes(RDO_mgL_9_adjusted, CTD_DO_9, colour=factor(time))) + 
  geom_point(alpha = 0.6) + 
  geom_abline(intercept = as.numeric(EXO2ctd_do_9_model$coefficients[1]),
              slope = as.numeric(EXO2ctd_do_9_model$coefficients[2]), colour="black") +
  geom_abline(intercept = 0, slope = 1, colour = "red")#1:1 line

#summary stats for chla model
print(EXO2ctd_do_9_model$coefficients[1])#intercept
print(EXO2ctd_do_9_model$coefficients[2])#slope
print(sd(EXO2ctd_do_9_model$residuals)) #standard deviation of the residuals
#THE REGRESSION MODEL AND ONE-ONE LINE ARE SO CLOSE, NO NEED TO CONVERT EXO --> CTD AT 9M!



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
#convert 1.6 m EXO DO to CTD DO

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





#convert 5 & 9 m DO sensor to CTD DO


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

#pull in FCR chem data from 2013-2020 from EDI
if(!file.exists("chemistry_HLW_edited.csv")){
  inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/199/8/da174082a3d924e989d3151924f9ef98" 
  infile1 <- paste0(getwd(),"/chemistry.csv")
  download.file(inUrl1,infile1,method="curl")
}

FCRchem <- readr::read_csv("chemistry_HLW_edited.csv") %>%
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
  dplyr::summarise(across(c("EXODO_mgL_1":"EXOfDOM_QSU_1"),mean)) 

#convert 1.6 EXO chla into CTD chla
fdom2doc <- dplyr::left_join(catwalk, FCRchem, by = c('time')) %>% 
  select(time, EXOfDOM_QSU_1, DOC_mgL) %>% 
  mutate(OGM_doc_total =)



#convert 1.6 m EXO DO to CTD DO
#convert 5 & 9 m DO sensor to CTD DO


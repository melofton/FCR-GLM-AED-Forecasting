setwd("../inputs")
wk_dir <- getwd()
remotes::install_github("FLARE-forecast/FLAREr")
remotes::install_github("FLARE-forecast/Rnoaa4cast")
library(FLAREr)


if(!file.exists("Met_final_2015_2020.csv")){
download.file("https://pasta.lternet.edu/package/data/eml/edi/389/5/3d1866fecfb8e17dc902c76436239431", 
              destfile = "Met_final_2015_2020.csv", 
              method="curl")
}

#File from Github
download.file("https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-metstation-data/FCRmet.csv", 
              destfile = "FCRmet.csv", 
              method="curl")


observed_met_file <- file.path(wk_dir, paste0("observed-met_fcre.nc"))

source("met_qaqc.R")
met_qaqc(realtime_file = "FCRmet.csv",
         qaqc_file = "Met_final_2015_2020.csv",
         cleaned_met_file_dir = wk_dir,
         input_file_tz = "EST",
         nldas = NULL)


#need to get the relevent variables from the config list and set them here
config <- list()
config$met$use_forecasted_met <- FALSE
config$run_config$start_datetime <- "2015-07-07"
config$run_config$end_datetime <- "2020-12-31"
config$run_config$forecast_start_datetime <- "2020-12-31"


met_out <- FLAREr::generate_glm_met_files(obs_met_file = observed_met_file,
                                          out_dir = wk_dir,
                                          forecast_dir = NULL,
                                          config = config)

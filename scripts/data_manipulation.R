library(dplyr)
library(readr)
library(readxl)
library(purrr)
library(stringr)
library(tidyr)

x <- list.dirs(path = "./rawData/SMEAR2/Larger") %>%
  lapply(function(x) {list.files(x, full.names = TRUE)}) %>%
  .[-1] %>%
  lapply(function(x) {lapply(x, read.csv) %>% bind_rows}) %>%
  lapply(function(x) {x$Time <- as.POSIXct(str_c(x$Year, "-", x$Month, "-", x$Day, " ", x$Hour, ":", x$Minute, ":00"), format="%Y-%m-%d %H:%M:%S", tz="UTC");x}) %>%
  #lapply(function(x) filter(x, Month > 6 & Month < 9)) %>%
  lapply(function(x) x[!names(x) %in% c("Minute","Second", "Hour", "Day", "Year", "Month")]) %>%
  lapply(function(x) arrange(x, Time)) %>%
  lapply(function(x){if(length(x$Time) > 464400) x[-c(1:4368),] else x}) %>%
  reduce(cbind) %>%
  .[unique(colnames(.))] 

path1 <- "./rawData/SMEAR2/Larger/smeardata_magicPAR.csv"
path2 <- "./rawData/SMEAR2/Larger/smeardata_magicDIFFUSEPAR.csv"

par <- read.csv(path1)
diffpar <- read.csv(path2)
par$HYY_META.PAR2_magic <-par$HYY_META.Glob * 2.06
par$HYY_META.diffPAR_magic <- diffpar$HYY_META.diffGLOB * 2.06
par$Time <- as.POSIXct(str_c(par$Year, "-", par$Month, "-", par$Day, " ", par$Hour, ":", par$Minute, ":00"), format="%Y-%m-%d %H:%M:%S", tz="UTC")
par <- par[!names(par) %in% c("Minute","Second", "Hour", "Day", "Year", "Month", "HYY_META.Glob")]

x <- merge(x, par, by="Time", all.x = TRUE)
x$HYY_META.PAR2 <- coalesce(x$HYY_META.PAR2, x$HYY_META.PAR2_magic)
x$HYY_META.diffPAR <- coalesce(x$HYY_META.diffPAR, x$HYY_META.diffPAR_magic)
x <- x[!names(x) %in% c("HYY_META.PAR2_magic", "HYY_META.diffPAR_magic")]

soil <- x %>% select(HYY_META.wsoil_B2, Time) %>%
  na.omit() %>%
  mutate(Time = format(Time, '%Y-%m-%d')) %>%
  aggregate(HYY_META.wsoil_B2 ~ Time, mean) %>%
  mutate(Time = as.POSIXct(format(as.POSIXct(Time), '%Y-%m-%d %H:%M:%S'))) %>%
  mutate(Time = lapply(Time, function(x) seq(x, by = "30 min", length = 48))) %>%
  unnest(cols = c(Time))

x <- merge(x, soil, by="Time", all.x = TRUE)
x$HYY_META.wsoil_B2 <- coalesce(x$HYY_META.wsoil_B2.x, x$HYY_META.wsoil_B2.y)
x <- x[!names(x) %in% c("HYY_META.wsoil_B2.x", "HYY_META.wsoil_B2.y")]

x <- na.omit(x[!names(x) %in% c("HYY_META.wsoil_C1", "HYY_META.wsoil_A", "HYY_META.wsoil_B1")]) %>%
  filter(HYY_EDDY233.Qc_gapf_NEE == 0) %>%
  filter(HYY_META.PAR2 > 10) %>%
  within(rm(HYY_EDDY233.Qc_gapf_NEE))

names(x)[names(x) == 'HYY_EDDY233.u_star'] <- 'FricVel'
names(x)[names(x) == 'HYY_META.T336'] <- 'AirTemp'
names(x)[names(x) == 'HYY_META.tsoil_B2'] <- 'SoilTempB'
names(x)[names(x) == 'HYY_META.tsoil_A'] <- 'SoilTempA'
names(x)[names(x) == 'HYY_META.RHTd'] <- 'RelHum'
names(x)[names(x) == 'HYY_META.PAR2'] <- 'PAR'
names(x)[names(x) == 'HYY_META.wsoil_B2'] <- 'SoilWatCont'
names(x)[names(x) == 'HYY_META.diffPAR'] <- 'DiffRad'
names(x)[names(x) == 'HYY_EDDY233.NEE'] <- 'NEE'

x["DiffuseFract"] <- x["DiffRad"] / x["PAR"]
x["DiffuseFract"] <- x['DiffuseFract']

e_s <- 611 * exp( (17.27 * x["AirTemp"]) / (237.3 + x["AirTemp"]) )
e_a <- e_s * (x["RelHum"]/100)
x["VaporPressureDeficit"] = e_s - e_a

write.csv(x, "./data/smear2008-18_allmonths.csv")



x <- list.files(path = "./rawData/SMEAR1/varrio", full.names = TRUE) %>%
  lapply(read.csv) %>%
  lapply(function(x) {x$Time <- as.POSIXct(str_c(x$Year, "-", x$Month, "-", x$Day, " ", x$Hour, ":", x$Minute, ":00"), format="%Y-%m-%d %H:%M:%S", tz="UTC");x}) %>%
  #lapply(function(x) filter(x, Month > 6 & Month < 9 & Year > 2012)) %>%
  lapply(function(x) x[!names(x) %in% c("Minute","Second", "Hour", "Day", "Year", "Month")]) %>%
  lapply(function(x) arrange(x, Time))  %>%
  reduce(merge)


# 
soil <- x %>% select(VAR_META.wsoil, Time) %>%
  na.omit() %>%
  mutate(Time = format(Time, '%Y-%m-%d')) %>%
  aggregate(VAR_META.wsoil ~ Time, mean) %>%
  mutate(Time = as.POSIXct(format(as.POSIXct(Time), '%Y-%m-%d %H:%M:%S'))) %>%
  mutate(Time = lapply(Time, function(x) seq(x, by = "30 min", length = 48))) %>%
  unnest(cols = c(Time))
# 
x <- merge(x, soil, by="Time", all.x = TRUE)
x$VAR_META.wsoil <- coalesce(x$VAR_META.wsoil.x, x$VAR_META.wsoil.y)
x <- x[!names(x) %in% c("VAR_META.wsoil.x", "VAR_META.wsoil.y")]

# x <- na.omit(x[!names(x) %in% c("")]) %>%
#   filter(VAR_META.PAR > 10) %>%
#   filter(VAR_EDDY.Qc_gapf_NEE == 0) %>%
#   within(rm(VAR_EDDY.Qc_gapf_NEE))
# 
# 
names(x)[names(x) == 'VAR_EDDY.u_star'] <- 'FricVel'
names(x)[names(x) == 'VAR_META.TDRY1'] <- 'AirTemp'
names(x)[names(x) == 'VAR_META.ST'] <- 'SoilTempA'
names(x)[names(x) == 'VAR_META.HUM_RH'] <- 'RelHum'
names(x)[names(x) == 'VAR_META.PAR'] <- 'PAR'
names(x)[names(x) == 'VAR_META.wsoil'] <- 'SoilWatCont'
names(x)[names(x) == 'VAR_META.diffPAR'] <- 'DiffRad'
names(x)[names(x) == 'VAR_EDDY.NEE'] <- 'NEE'

x["DiffuseFract"] <- x["DiffRad"] / x["PAR"]
x["DiffuseFract"] <- x['DiffuseFract']

e_s <- 611 * exp( (17.27 * x["AirTemp"]) / (237.3 + x["AirTemp"]) )
e_a <- e_s * (x["RelHum"]/100)
x["VaporPressureDeficit"] = e_s - e_a


write.csv(x, "./data/varrio_2015-2019_allmonths.csv")

path = "./rawData/SMEAR2/hyytiala2019+/"
x <- list.files(path, full.names = TRUE)[-1] %>%
  lapply(read.csv) %>%
  lapply(function(x) {x$Time <- as.POSIXct(str_c(x$Year, "-", x$Month, "-", x$Day, " ", x$Hour, ":", x$Minute, ":00"), format="%Y-%m-%d %H:%M:%S", tz="UTC");x}) %>%
  #lapply(function(x) filter(x, Month > 6 & Month < 9)) %>%
  lapply(function(x) x[!names(x) %in% c("Minute","Second", "Hour", "Day", "Year", "Month")]) %>%
  lapply(function(x) arrange(x, Time))  %>%
  reduce(cbind) %>%
  .[unique(colnames(.))]

soil <- x %>% select(HYY_META.wsoil_B2, Time) %>%
  na.omit() %>%
  mutate(Time = format(Time, '%Y-%m-%d')) %>%
  aggregate(HYY_META.wsoil_B2 ~ Time, mean) %>%
  mutate(Time = as.POSIXct(format(as.POSIXct(Time), '%Y-%m-%d %H:%M:%S'))) %>%
  mutate(Time = lapply(Time, function(x) seq(x, by = "30 min", length = 48))) %>%
  unnest(cols = c(Time))

x <- merge(x, soil, by="Time", all.x = TRUE)
x$HYY_META.wsoil_B2 <- coalesce(x$HYY_META.wsoil_B2.x, x$HYY_META.wsoil_B2.y)
x <- x[!names(x) %in% c("HYY_META.wsoil_B2.x", "HYY_META.wsoil_B2.y")]

x <- na.omit(x[!names(x) %in% c("HYY_META.tsoil_B2")]) %>%
  filter(HYY_EDDY233.Qc_gapf_NEE == 0) %>%
  filter(HYY_META.PAR2 > 10) %>%
  within(rm(HYY_EDDY233.Qc_gapf_NEE))

names(x)[names(x) == 'HYY_EDDYMAST.u_star_270'] <- 'FricVel'
names(x)[names(x) == 'HYY_META.T336'] <- 'AirTemp'
names(x)[names(x) == 'HYY_META.tsoil_B2'] <- 'SoilTempB'
names(x)[names(x) == 'HYY_META.tsoil_A'] <- 'SoilTempA'
names(x)[names(x) == 'HYY_META.RHTd'] <- 'RelHum'
names(x)[names(x) == 'HYY_META.PAR2'] <- 'PAR'
names(x)[names(x) == 'HYY_META.wsoil_B2'] <- 'SoilWatCont'
names(x)[names(x) == 'HYY_META.diffPAR'] <- 'DiffRad'
names(x)[names(x) == 'HYY_EDDY233.NEE'] <- 'NEE'

x["DiffuseFract"] <- x["DiffRad"] / x["PAR"]
x["DiffuseFract"] <- x['DiffuseFract']

e_s <- 611 * exp((17.27 * x["AirTemp"]) / (237.3 + x["AirTemp"]))
e_a <- e_s * (x["RelHum"]/100)
x["VaporPressureDeficit"] = e_s - e_a

write.csv(x, "./data/smear2019-21_allmonths.csv")

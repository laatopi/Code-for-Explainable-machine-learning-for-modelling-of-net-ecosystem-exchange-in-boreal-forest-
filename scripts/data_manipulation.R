# Load required libraries
library(dplyr)
library(readr)
library(readxl)
library(purrr)
library(stringr)
library(tidyr)

# Load configuration file for paths
source("config.R")

# SMEAR2 Processing
smear2_path <- file.path(raw_data_path, "SMEAR2/Larger/")
x <- list.dirs(path = smear2_path, recursive = FALSE) %>%
  lapply(function(dir) list.files(dir, full.names = TRUE)) %>%
  .[-1] %>%
  lapply(function(files) lapply(files, read.csv) %>% bind_rows()) %>%
  lapply(function(df) {
    df$Time <- as.POSIXct(str_c(df$Year, "-", df$Month, "-", df$Day, " ", df$Hour, ":", df$Minute, ":00"), format="%Y-%m-%d %H:%M:%S", tz="UTC")
    df
  }) %>%
  lapply(function(df) select(df, -Minute, -Second, -Hour, -Day, -Year, -Month)) %>%
  lapply(function(df) arrange(df, Time)) %>%
  lapply(function(df) if (nrow(df) > 464400) df[-c(1:4368),] else df) %>%
  reduce(cbind) %>%
  .[unique(colnames(.))]

# Load additional PAR data
path1 <- file.path(smear2_path, "smeardata_magicPAR.csv")
path2 <- file.path(smear2_path, "smeardata_magicDIFFUSEPAR.csv")

par <- read.csv(path1)
diffpar <- read.csv(path2)
par <- par %>%
  mutate(
    HYY_META.PAR2_magic = HYY_META.Glob * 2.06,
    HYY_META.diffPAR_magic = diffpar$HYY_META.diffGLOB * 2.06,
    Time = as.POSIXct(str_c(Year, "-", Month, "-", Day, " ", Hour, ":", Minute, ":00"), format="%Y-%m-%d %H:%M:%S", tz="UTC")
  ) %>%
  select(-Minute, -Second, -Hour, -Day, -Year, -Month, -HYY_META.Glob)

# Merge PAR data into main data
x <- merge(x, par, by = "Time", all.x = TRUE) %>%
  mutate(
    HYY_META.PAR2 = coalesce(HYY_META.PAR2, HYY_META.PAR2_magic),
    HYY_META.diffPAR = coalesce(HYY_META.diffPAR, HYY_META.diffPAR_magic)
  ) %>%
  select(-HYY_META.PAR2_magic, -HYY_META.diffPAR_magic)

# Soil moisture processing for SMEAR2
soil <- x %>%
  select(HYY_META.wsoil_B2, Time) %>%
  na.omit() %>%
  mutate(Time = format(Time, '%Y-%m-%d')) %>%
  aggregate(HYY_META.wsoil_B2 ~ Time, mean) %>%
  mutate(Time = as.POSIXct(Time, format = '%Y-%m-%d %H:%M:%S')) %>%
  unnest(cols = list(Time = seq(Time, by = "30 min", length.out = 48)))

x <- merge(x, soil, by = "Time", all.x = TRUE) %>%
  mutate(HYY_META.wsoil_B2 = coalesce(HYY_META.wsoil_B2.x, HYY_META.wsoil_B2.y)) %>%
  select(-HYY_META.wsoil_B2.x, -HYY_META.wsoil_B2.y)

# Filter and rename columns for SMEAR2
x <- x %>%
  filter(HYY_EDDY233.Qc_gapf_NEE == 0, HYY_META.PAR2 > 10) %>%
  select(-HYY_EDDY233.Qc_gapf_NEE, -HYY_META.wsoil_C1, -HYY_META.wsoil_A, -HYY_META.wsoil_B1) %>%
  rename(
    FricVel = HYY_EDDY233.u_star,
    AirTemp = HYY_META.T336,
    SoilTempB = HYY_META.tsoil_B2,
    SoilTempA = HYY_META.tsoil_A,
    RelHum = HYY_META.RHTd,
    PAR = HYY_META.PAR2,
    SoilWatCont = HYY_META.wsoil_B2,
    DiffRad = HYY_META.diffPAR,
    NEE = HYY_EDDY233.NEE
  ) %>%
  mutate(
    DiffuseFract = DiffRad / PAR,
    VaporPressureDeficit = 611 * exp((17.27 * AirTemp) / (237.3 + AirTemp)) - 
      (611 * exp((17.27 * AirTemp) / (237.3 + AirTemp)) * (RelHum / 100))
  )

# Save processed SMEAR2 data
write.csv(x, file.path(base_output_path, "smear2008-18_allmonths.csv"))

# Process SMEAR1 Data
smear1_path <- file.path(raw_data_path, "SMEAR1/varrio/")
x <- list.files(path = smear1_path, full.names = TRUE) %>%
  lapply(read.csv) %>%
  lapply(function(df) {
    df$Time <- as.POSIXct(str_c(df$Year, "-", df$Month, "-", df$Day, " ", df$Hour, ":", df$Minute, ":00"), format="%Y-%m-%d %H:%M:%S", tz="UTC")
    df
  }) %>%
  lapply(function(df) select(df, -Minute, -Second, -Hour, -Day, -Year, -Month)) %>%
  lapply(function(df) arrange(df, Time)) %>%
  reduce(merge)

# Soil processing for SMEAR1
soil <- x %>%
  select(VAR_META.wsoil, Time) %>%
  na.omit() %>%
  mutate(Time = format(Time, '%Y-%m-%d')) %>%
  aggregate(VAR_META.wsoil ~ Time, mean) %>%
  mutate(Time = as.POSIXct(Time, format = '%Y-%m-%d %H:%M:%S')) %>%
  unnest(cols = list(Time = seq(Time, by = "30 min", length.out = 48)))

x <- merge(x, soil, by = "Time", all.x = TRUE) %>%
  mutate(VAR_META.wsoil = coalesce(VAR_META.wsoil.x, VAR_META.wsoil.y)) %>%
  select(-VAR_META.wsoil.x, -VAR_META.wsoil.y) %>%
  rename(
    FricVel = VAR_EDDY.u_star,
    AirTemp = VAR_META.TDRY1,
    SoilTempA = VAR_META.ST,
    RelHum = VAR_META.HUM_RH,
    PAR = VAR_META.PAR,
    SoilWatCont = VAR_META.wsoil,
    DiffRad = VAR_META.diffPAR,
    NEE = VAR_EDDY.NEE
  ) %>%
  mutate(
    DiffuseFract = DiffRad / PAR,
    VaporPressureDeficit = 611 * exp((17.27 * AirTemp) / (237.3 + AirTemp)) - 
      (611 * exp((17.27 * AirTemp) / (237.3 + AirTemp)) * (RelHum / 100))
  )

# Save processed SMEAR1 data
write.csv(x, file.path(base_output_path, "varrio_2015-2019_allmonths.csv"))

# Process Post-2019 SMEAR2 Data
smear2_2019_path <- file.path(raw_data_path, "SMEAR2/hyytiala2019+/")
x <- list.files(smear2_2019_path, full.names = TRUE)[-1] %>%
  lapply(read.csv) %>%
  lapply(function(df) {
    df$Time <- as.POSIXct(str_c(df$Year, "-", df$Month, "-", df$Day, " ", df$Hour, ":", df$Minute, ":00"), format="%Y-%m-%d %H:%M:%S", tz="UTC")
    df
  }) %>%
  lapply(function(df) select(df, -Minute, -Second, -Hour, -Day, -Year, -Month)) %>%
  lapply(function(df) arrange(df, Time)) %>%
  reduce(cbind) %>%
  .[unique(colnames(.))]

# Soil moisture processing for Post-2019 SMEAR2 data
soil <- x %>%
  select(HYY_META.wsoil_B2, Time) %>%
  na.omit() %>%
  mutate(Time = format(Time, '%Y-%m-%d')) %>%
  aggregate(HYY_META.wsoil_B2 ~ Time, mean) %>%
  mutate(Time = as.POSIXct(Time, format = '%Y-%m-%d %H:%M:%S')) %>%
  unnest(cols = list(Time = seq(Time, by = "30 min", length.out = 48)))

x <- merge(x, soil, by = "Time", all.x = TRUE) %>%
  mutate(HYY_META.wsoil_B2 = coalesce(HYY_META.wsoil_B2.x, HYY_META.wsoil_B2.y)) %>%
  select(-HYY_META.wsoil_B2.x, -HYY_META.wsoil_B2.y) %>%
  filter(HYY_EDDY233.Qc_gapf_NEE == 0, HYY_META.PAR2 > 10) %>%
  rename(
    FricVel = HYY_EDDYMAST.u_star_270,
    AirTemp = HYY_META.T336,
    SoilTempB = HYY_META.tsoil_B2,
    SoilTempA = HYY_META.tsoil_A,
    RelHum = HYY_META.RHTd,
    PAR = HYY_META.PAR2,
    SoilWatCont = HYY_META.wsoil_B2,
    DiffRad = HYY_META.diffPAR,
    NEE = HYY_EDDY233.NEE
  ) %>%
  mutate(
    DiffuseFract = DiffRad / PAR,
    VaporPressureDeficit = 611 * exp((17.27 * AirTemp) / (237.3 + AirTemp)) - 
      (611 * exp((17.27 * AirTemp) / (237.3 + AirTemp)) * (RelHum / 100))
  )

# Save Post-2019 SMEAR2 processed data
write.csv(x, file.path(base_output_path, "smear2019-21_allmonths.csv"))

# 10_weather_data_cleaning.R ----
# Downloads DWD HYRAS daily temperature (TAS) and precipitation (PR) rasters for 2024,
# extracts Apr-Jul daily values at the 20 Whinchat breeding sites (bilinear interpolation)
# and aggregates them to Jun-Jul breeding summaries
# Site coordinates are projected (ETRS89/UTM32N, EPSG:25832), the HYRAS grid is ETRS89-LAEA (EPSG:3035)
# TAS is a calendar day mean, PR is a daily sum from 06:00 UTC to 06:00 UTC of the following day
# Requires: data_clean/site_ref.rds (produced by script 00)

# ---- packages ----
library("terra")
library("sf")
library("tidyverse")
library("readxl")
library("writexl")
library("here")

# ---- paths ----
sites_data_raw <- here("data_raw/mastertable.xlsx")
hyras_dir <- here("data_raw/hyras")
output_dir <- here("data_clean")
dir.create(hyras_dir, showWarnings = FALSE, recursive = TRUE)

# ---- settings ----
study_year <- 2024L
extract_months <- 4:7
breeding_months <- 6:7
dry_day_mm <- 1
heavy_rain_mm <- 15
options(timeout = 1800)

# ---- site reference ----
site_ref <- readRDS(here("data_clean/site_ref.rds"))
breeding_ids <- site_ref |> filter(breeding == "ja") |> pull(Sample_ID)

# ---- sites ----
sites <- read_excel(sites_data_raw) |>
  filter(breeding == "ja") |>
  left_join(select(site_ref, ID, Sample_ID), by = "ID", relationship = "one-to-one") |>
  rename(Site_ID = Sample_ID) |>
  select(Site_ID, Lon, Lat)

stopifnot(
  nrow(sites) == 20,
  !anyDuplicated(sites$Site_ID),
  setequal(sites$Site_ID, breeding_ids),
  !anyNA(sites[c("Lon", "Lat")]),
  between(sites$Lon, 4e5, 6e5), # UTM32N easting range for Hesse
  between(sites$Lat, 5.4e6, 5.7e6) # UTM32N northing range for Hesse
)

sites <- st_as_sf(sites, coords = c("Lon", "Lat"), crs = 25832)

# ---- HYRAS download helpers ----
hyras_url <- function(var, yr) {
  base <- "https://opendata.dwd.de/climate_environment/CDC/grids_germany/daily/hyras_de/"
  sub <- if (var == "tas") "air_temperature_mean" else "precipitation"
  sprintf("%s%s/%s_hyras_1_%d_v6-0_de.nc", base, sub, var, yr)
}

# download to a temporary file first so an interrupted download is never mistaken for a complete file
download_hyras <- function(var, yr) {
  dest <- file.path(hyras_dir, sprintf("%s_%d.nc", var, yr))
  if (!file.exists(dest)) {
    part <- paste0(dest, ".part")
    status <- download.file(hyras_url(var, yr), part, mode = "wb", quiet = TRUE)
    stopifnot(status == 0, file.size(part) > 0)
    file.rename(part, dest)
  }
  dest
}

# ---- raw rasters ----
r_tas <- rast(download_hyras("tas", study_year))
r_pr <- rast(download_hyras("pr", study_year))

stopifnot(
  compareGeom(r_tas, r_pr, stopOnError = FALSE),
  identical(as.Date(time(r_tas)), as.Date(time(r_pr))),
  all(year(as.Date(time(r_tas))) == study_year)
)

# ---- extract Apr-Jul daily values ----
dates <- as.Date(time(r_tas))
sel <- month(dates) %in% extract_months
r_tas <- r_tas[[sel]]
r_pr <- r_pr[[sel]]
dates <- dates[sel]

pts <- project(vect(sites), crs(r_tas))
v_tas <- as.matrix(terra::extract(r_tas, pts, method = "bilinear")[, -1])
v_pr <- as.matrix(terra::extract(r_pr, pts, method = "bilinear")[, -1])

# sites outside the raster would return NA, so this must be checked before aggregating
stopifnot(!anyNA(v_tas), !anyNA(v_pr))

weather_daily <- map(seq_len(nrow(sites)), \(i) tibble(
  Site_ID = sites$Site_ID[i],
  Date = dates,
  temp_mean_c = v_tas[i, ],
  precip_mm = v_pr[i, ]
)) |> list_rbind() |>
  arrange(Site_ID, Date)

# ---- check (daily) ----
stopifnot(
  nrow(weather_daily) == nrow(sites) * length(dates),
  setequal(weather_daily$Site_ID, breeding_ids),
  !anyNA(weather_daily),
  !anyDuplicated(weather_daily[c("Site_ID", "Date")]),
  between(weather_daily$temp_mean_c, -30, 45),
  all(weather_daily$precip_mm >= 0 & weather_daily$precip_mm < 400),
  all(count(filter(weather_daily, month(Date) %in% breeding_months), Site_ID)$n == 61)
)

# ---- breeding summary: Jun-Jul ----
weather_breeding <- weather_daily |>
  filter(month(Date) %in% breeding_months) |>
  group_by(Site_ID) |>
  summarise(
    temp_mean_breeding = mean(temp_mean_c),
    max_daily_precip_breeding = max(precip_mm),
    n_drydays_breeding = sum(precip_mm < dry_day_mm),
    n_heavy_rain_days_breeding = sum(precip_mm >= heavy_rain_mm),
    .groups = "drop"
  )

# ---- check (breeding summary) ----
stopifnot(
  nrow(weather_breeding) == nrow(sites),
  setequal(weather_breeding$Site_ID, breeding_ids),
  !anyDuplicated(weather_breeding$Site_ID),
  !anyNA(weather_breeding),
  all(weather_breeding$n_drydays_breeding <= 61),
  all(weather_breeding$n_heavy_rain_days_breeding <= 61),
  all((weather_breeding$n_heavy_rain_days_breeding > 0) ==
        (weather_breeding$max_daily_precip_breeding >= heavy_rain_mm))
)

# ---- save ----
saveRDS(weather_daily, file.path(output_dir, "weather_daily.rds"))
saveRDS(weather_breeding, file.path(output_dir, "weather_breeding.rds"))

 write_xlsx(weather_daily, file.path(output_dir, "weather_daily.xlsx"))
 write_xlsx(weather_breeding, file.path(output_dir, "weather_breeding.xlsx"))
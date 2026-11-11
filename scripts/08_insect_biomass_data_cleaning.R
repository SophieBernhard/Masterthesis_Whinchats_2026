# 08_insect_biomass_data_cleaning.R ----
# Load, clean, and aggregate insect biomass data per site and phase
# Arrival = 2 trap events (17 and 31 May), breeding = 4 trap events (14 Jun, 1 Jul, 16 Jul, 29 Jul)
# Weighings of damaged traps are set to NA and means are taken over valid weighings only
# Requires: data_clean/site_ref.rds (produced by script 00)

# ---- packages ----
library("readxl")
library("tidyverse")
library("writexl")
library("here")

# ---- paths ----
biomass_data_raw <- here("data_raw/insect-biomass_land-use_weather_masterthesi_Schmitt.xlsx")
output_dir <- here("data_clean")

# ---- site reference ----
site_ref <- readRDS(here("data_clean/site_ref.rds"))

# ---- settings ----
# Phase assignment by trap date (day.month)
phase_map <- c(
  "17.5" = "arrival",
  "31.5" = "arrival",
  "14.6" = "breeding",
  "1.7" = "breeding",
  "16.7" = "breeding",
  "29.7" = "breeding"
)

# ---- raw data ----
biomass_raw <- read_excel(biomass_data_raw) |>
  rename_with(~ "weight_g", .cols = matches("^Weight", ignore.case = TRUE)) |>
  rename_with(~ "date_num", .cols = matches("^Date", ignore.case = TRUE)) |>
  rename_with(~ "trap_damage", .cols = matches("damage", ignore.case = TRUE)) |>
  select(Site_ID = Trap, weight_g, date_num, trap_damage) |>
  filter(Site_ID %in% site_ref$Sample_ID) |>
  mutate(
    weight_g = as.numeric(weight_g),
    trap_damage = tolower(trimws(trap_damage)),
    date_chr = as.character(date_num)
  )

# ---- check (raw data) ----
# Every site has all six trap events, all dates are known, no weighing is missing
stopifnot(
  setequal(biomass_raw$Site_ID, site_ref$Sample_ID),
  nrow(biomass_raw) == length(phase_map) * n_distinct(site_ref$Sample_ID),
  !anyDuplicated(biomass_raw[c("Site_ID", "date_chr")]),
  all(biomass_raw$date_chr %in% names(phase_map)),
  all(biomass_raw$trap_damage %in% c("yes", "no")),
  !anyNA(biomass_raw$weight_g),
  all(biomass_raw$weight_g >= 0)
)

# ---- aggregate ----
biomass <- biomass_raw |>
  mutate(
    weight_g = if_else(trap_damage == "yes", NA_real_, weight_g),
    phase = phase_map[date_chr]
  ) |>
  group_by(Site_ID, phase) |>
  summarise(
    mean_biomass = mean(weight_g, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(mean_biomass = if_else(is.nan(mean_biomass), NA_real_, mean_biomass)) |>
  pivot_wider(
    names_from = phase,
    values_from = mean_biomass,
    names_prefix = "biomass_mean_"
  ) |>
  arrange(Site_ID)

# ---- check (result) ----
# LDK17 is the only site without a valid arrival weighing (both arrival traps damaged)
stopifnot(
  nrow(biomass) == n_distinct(site_ref$Sample_ID),
  setequal(biomass$Site_ID, site_ref$Sample_ID),
  !anyDuplicated(biomass$Site_ID),
  setequal(names(biomass), c("Site_ID", "biomass_mean_arrival", "biomass_mean_breeding")),
  setequal(biomass$Site_ID[is.na(biomass$biomass_mean_arrival)], "LDK17"),
  !anyNA(biomass$biomass_mean_breeding),
  all(biomass$biomass_mean_arrival > 0, na.rm = TRUE),
  all(biomass$biomass_mean_breeding > 0)
)

# ---- save ----
saveRDS(biomass, file.path(output_dir, "biomass.rds"))
 write_xlsx(biomass, file.path(output_dir, "biomass.xlsx"))
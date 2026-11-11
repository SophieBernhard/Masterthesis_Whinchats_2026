# 09_land_use_percentage_data_cleaning.R ----
# Load and clean land-use composition (% of area in a 100 m radius) per site (n = 40 plots)
# Requires: data_clean/site_ref.rds (produced by script 00)

# ---- packages ----
library("readxl")
library("tidyverse")
library("writexl")
library("here")

# ---- paths ----
landuse_data_raw <- here("data_raw/insect-biomass_land-use_weather_masterthesi_Schmitt.xlsx")
output_dir <- here("data_clean")

# ---- site reference ----
site_ref <- readRDS(here("data_clean/site_ref.rds"))

# ---- raw records ----
# six trap events per site, land-use columns are constant within a site
landuse_raw <- read_excel(landuse_data_raw) |>
  rename_with(~ str_remove_all(.x, "[\r\n\\[\\]%]") |> str_trim()) |>
  select(Site_ID = Trap, Extensive, Intensive, Arableland, Seminat, Settlement, Trees)

stopifnot(
  nrow(landuse_raw) == 6 * n_distinct(site_ref$Sample_ID),
  setequal(landuse_raw$Site_ID, site_ref$Sample_ID),
  all(count(landuse_raw, Site_ID)$n == 6),
  all(count(distinct(landuse_raw), Site_ID)$n == 1)
)

# ---- site level table ----
landuse_site <- landuse_raw |>
  distinct() |>
  arrange(Site_ID)

# ---- check ----
# all six classes are percentages and sum to 100 within rounding error
class_cols <- c("Extensive", "Intensive", "Arableland", "Seminat", "Settlement", "Trees")
stopifnot(
  nrow(landuse_site) == 40,
  !anyDuplicated(landuse_site$Site_ID),
  !anyNA(landuse_site),
  all(landuse_site[class_cols] >= 0 & landuse_site[class_cols] <= 100),
  all(abs(rowSums(landuse_site[class_cols]) - 100) < 0.01)
)

# ---- analysis table ----
landuse_clean <- landuse_site |>
  select(
    Site_ID,
    pct_extensive = Extensive,
    pct_intensive = Intensive,
    pct_arable = Arableland,
    pct_seminat = Seminat,
    pct_settlement = Settlement
  )

stopifnot(
  nrow(landuse_clean) == 40,
  setequal(names(landuse_clean), c("Site_ID", "pct_extensive", "pct_intensive",
                                   "pct_arable", "pct_seminat", "pct_settlement"))
)

# ---- save ----
saveRDS(landuse_clean, file.path(output_dir, "land_use.rds"))
 write_xlsx(landuse_clean, file.path(output_dir, "land_use.xlsx"))
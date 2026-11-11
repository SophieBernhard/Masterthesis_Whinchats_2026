# 04_perches_breeding_data_cleaning.R ----
# Clean perch survey data from the breeding phase (June and July) and convert
# text-coded abundance categories to interval midpoints
# Requires: data_clean/site_ref.rds (produced by script 00)

# ---- packages ----
library("readxl")
library("tidyverse")
library("writexl")
library("here")

# ---- paths ----
perch_data_raw <- here("data_raw/04_Anlage C_Daten_Gesamt.xlsx")
output_dir <- here("data_clean")

# ---- site reference ----
site_ref <- readRDS(here("data_clean/site_ref.rds"))

# ---- midpoint lookup ----
# Absolute counts are kept, ranges are replaced by their midpoint, "> 50" is set
# to the midpoint of "50-100" and the open top category "> 500" is fixed at 600
perch_lookup <- c(
  "0" = 0L, "10" = 10L, "20" = 20L, "20-50" = 35L, "50-100" = 75L,
  ">50" = 75L, "100-200" = 150L, "200-500" = 350L, ">500" = 600L
)

# ---- read raw data ----
perches_raw <- read_excel(
  perch_data_raw,
  sheet = "B-1-Sitzwarten-Brutzeit",
  .name_repair = "unique_quiet"
) |>
  select(PlotID, raw_category = `Anzahl potenzieller Sitzwarten (Vegetation)`) |>
  mutate(
    PlotID_clean = str_remove_all(PlotID, " "),
    raw_category = str_remove_all(as.character(raw_category), " ")
  )

stopifnot(
  all(perches_raw$raw_category %in% names(perch_lookup)),
  !anyDuplicated(perches_raw$PlotID_clean),
  all(perches_raw$PlotID_clean %in% site_ref$PlotID_clean)
)

# ---- analysis table ----
perches_breeding <- perches_raw |>
  left_join(
    select(site_ref, Sample_ID, PlotID_clean),
    by = "PlotID_clean",
    relationship = "many-to-one"
  ) |>
  mutate(n_perches_mid = unname(perch_lookup[raw_category])) |>
  select(Site_ID = Sample_ID, raw_category, n_perches_mid) |>
  arrange(Site_ID)

stopifnot(
  nrow(perches_breeding) == 40,
  !anyDuplicated(perches_breeding$Site_ID),
  !anyNA(perches_breeding)
)

# ---- save ----
saveRDS(perches_breeding, file.path(output_dir, "perches_breeding.rds"))
 write_xlsx(perches_breeding, file.path(output_dir, "perches_breeding.xlsx"))
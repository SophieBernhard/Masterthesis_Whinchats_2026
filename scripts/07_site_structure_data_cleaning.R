# 07_site_structure_data_cleaning.R ----
# Clean and select habitat structure variables from drone LiDAR data (n = 40 plots)
# Drone data only (equivalent to flight data, temporally more appropriate; Steen, 2025)
# Raw values are saved, transformations are applied in script 13
# Requires: data_clean/site_ref.rds (produced by script 00)

# ---- packages ----
library("readxl")
library("tidyverse")
library("writexl")
library("here")

# ---- paths ----
struct_data_raw <- here("data_raw/mastertable.xlsx")
output_dir <- here("data_clean")

# ---- site reference ----
site_ref <- readRDS(here("data_clean/site_ref.rds"))

# ---- raw records ----
# one row per plot, German raw names are kept in the file
# kraut_sd = SD of herb layer height, sitzwarten_25sqm = perches per 25 m2
struct_raw <- read_excel(struct_data_raw, .name_repair = "unique_quiet")

stopifnot(
  nrow(struct_raw) == 40,
  !anyDuplicated(struct_raw$ID),
  !anyNA(struct_raw$ID),
  setequal(struct_raw$ID, site_ref$ID)
)

# ---- join site IDs and select ----
site_structure <- struct_raw |>
  left_join(
    site_ref |> select(ID, Sample_ID, breeding_ref = breeding),
    by = "ID",
    relationship = "one-to-one"
  )

# raw breeding status must match the site reference
stopifnot(
  !anyNA(site_structure$Sample_ID),
  all(site_structure$breeding == site_structure$breeding_ref)
)

site_structure <- site_structure |>
  select(
    Site_ID = Sample_ID,
    viewshed = viewshed_drone,
    herb_h_sd = kraut_sd,
    mean_h = mean_h_drone,
    max_h = max_h_drone,
    n_perches_25sqm = sitzwarten_25sqm,
    terrain_sd = sd_dgm_drone,
    rumple = rumple_drone
  ) |>
  arrange(Site_ID)

# ---- check ----
stopifnot(
  nrow(site_structure) == 40,
  !anyDuplicated(site_structure$Site_ID),
  setequal(site_structure$Site_ID, site_ref$Sample_ID),
  !anyNA(site_structure),
  all(site_structure$viewshed >= 0 & site_structure$viewshed <= 1),
  all(site_structure$herb_h_sd >= 0),
  all(site_structure$mean_h > 0),
  all(site_structure$max_h >= site_structure$mean_h),
  all(site_structure$n_perches_25sqm >= 0),
  all(site_structure$n_perches_25sqm == round(site_structure$n_perches_25sqm)),
  all(site_structure$terrain_sd >= 0),
  all(site_structure$rumple >= 1)
)

site_structure <- site_structure |>
  mutate(n_perches_25sqm = as.integer(n_perches_25sqm))

# ---- save ----
saveRDS(site_structure, file.path(output_dir, "site_structure.rds"))
 write_xlsx(site_structure, file.path(output_dir, "site_structure.xlsx"))
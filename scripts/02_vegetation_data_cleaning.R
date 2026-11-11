# 02_vegetation_data_cleaning.R ----
# Clean raw Braun-Blanquet vegetation records and compute cover-weighted
# diversity metrics per site
# Requires: data_clean/site_ref.rds  (produced by script 00)

# ---- packages ----
library("readxl")
library("tidyverse")
library("writexl")
library("here")

# ---- paths ----
veg_data_raw <- here("data_raw/04_Anlage C_Daten_Gesamt.xlsx")
output_dir <- here("data_clean")

# ---- BB-scale lookup (class cover %) ----
# Turboveg defaults for the extended scale, see goeveg::scale_tabs$braun.blanquet2
# the single "-" record is treated as "r"
bb_lookup <- c("-" = 1, "r" = 1, "+" = 2, "1" = 3, "2m" = 4,
               "2a" = 8, "2b" = 18, "3" = 38, "4" = 63, "5" = 88)

# ---- site reference ----
site_ref <- readRDS(here("data_clean/site_ref.rds"))

# ---- raw vegetation data ----
veg_raw <- read_excel(
  veg_data_raw,
  sheet = "A-2-Vegetationsaufnahme",
  .name_repair = "unique_quiet"
) |>
  select(PlotID_clean = "PlotID", Species = "Art", BB_Code = "Deckungsgrad*") |>
  mutate(
    PlotID_clean = str_remove_all(PlotID_clean, " "),
    Species = str_squish(Species),
    BB_Code = as.character(BB_Code)
  )

stopifnot(
  all(veg_raw$BB_Code %in% names(bb_lookup)),
  all(veg_raw$PlotID_clean %in% site_ref$PlotID_clean)
)

# ---- join site IDs, convert BB codes, resolve duplicates ----
# duplicate plot-species records keep the highest cover class
vegetation <- veg_raw |>
  left_join(
    select(site_ref, Sample_ID, PlotID_clean),
    by = "PlotID_clean",
    relationship = "many-to-one"
  ) |>
  mutate(Coverage = bb_lookup[BB_Code]) |>
  group_by(Sample_ID, Species) |>
  slice_max(Coverage, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(Site_ID = Sample_ID, Species, BB_Code, Coverage) |>
  arrange(Site_ID, Species)

# ---- cover-weighted diversity per site ----
veg_diversity <- vegetation |>
  group_by(Site_ID) |>
  summarise(
    veg_richness = n_distinct(Species),
    veg_shannon = -sum((Coverage / sum(Coverage)) * log(Coverage / sum(Coverage))),
    .groups = "drop"
  ) |>
  mutate(veg_evenness = veg_shannon / log(veg_richness))

stopifnot(
  nrow(veg_diversity) == 40,
  !anyNA(veg_diversity),
  all(veg_diversity$veg_richness > 1)
)

# ---- save ----
saveRDS(vegetation, file.path(output_dir, "vegetation.rds"))
saveRDS(veg_diversity, file.path(output_dir, "veg_diversity.rds"))

 write_xlsx(vegetation,    file.path(output_dir, "vegetation.xlsx"))
 write_xlsx(veg_diversity, file.path(output_dir, "veg_diversity.xlsx"))
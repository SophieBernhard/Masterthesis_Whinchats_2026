# 03_perches_may_data_cleaning.R ----
# Clean perch survey data from the arrival phase (May) and keep the absolute
# count of potential perches per site
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

# ---- raw perch counts ----
perches_raw <- read_excel(
  perch_data_raw,
  sheet = "B-1-Sitzwarten-Ankunft",
  .name_repair = "unique_quiet"
) |>
  select(PlotID, n_perches = `Anzahl potenzieller Sitzwarten`) |>
  mutate(
    PlotID_clean = str_remove_all(PlotID, " "),
    n_perches = as.numeric(n_perches)
  )

# ---- correction ----
# VB21K is recorded with 0 perches, but a field note documents one tree at the
# site margin, so the count is set to 1
vb21k <- which(perches_raw$PlotID_clean == "VB21K")
stopifnot(
  length(vb21k) == 1,
  perches_raw$n_perches[vb21k] == 0,
  "VB21K" %in% site_ref$PlotID_clean
)
perches_raw$n_perches[vb21k] <- 1

# ---- join site IDs and check ----
perches_may <- perches_raw |>
  left_join(
    select(site_ref, Sample_ID, PlotID_clean),
    by = "PlotID_clean",
    relationship = "many-to-one"
  ) |>
  filter(!is.na(Sample_ID)) |>
  select(Site_ID = Sample_ID, n_perches) |>
  arrange(Site_ID)

stopifnot(
  nrow(perches_may) == 40,
  !anyDuplicated(perches_may$Site_ID),
  !anyNA(perches_may),
  all(perches_may$n_perches >= 0),
  all(perches_may$n_perches == round(perches_may$n_perches))
)
perches_may <- perches_may |> mutate(n_perches = as.integer(n_perches))

# ---- save ----
saveRDS(perches_may, file.path(output_dir, "perches_may.rds"))
 write_xlsx(perches_may, file.path(output_dir, "perches_may.xlsx"))

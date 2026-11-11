# 05_pesticides_data_cleaning.R ----
# Clean raw pesticide detection data and compute season-wide detection counts
# per breeding site (a compound counts as present if it was detected in any
# sample type on any survey date)
# Requires: data_clean/site_ref.rds (produced by script 00)

# ---- packages ----
library("readxl")
library("tidyverse")
library("writexl")
library("here")

# ---- paths ----
pest_data_raw <- here("data_raw/Braunkehlchen Total number of pesticides by month.xlsx")
output_dir <- here("data_clean")

# ---- settings ----
# EtOH = direct ethanol, eEtOH = evaporated ethanol, PG = propylene glycol
# "EtOH e" is a typo for "eEtOH" in the raw file (LDK15, 17 May)
sample_types <- c("EtOH" = "EtOH", "eEtOH" = "eEtOH", "EtOH e" = "eEtOH", "PG" = "PG")

# ---- site reference ----
site_ref <- readRDS(here("data_clean/site_ref.rds"))
breeding_ids <- site_ref |> filter(breeding == "ja") |> pull(Sample_ID)

# ---- pesticide names and category groups (rows 1 and 2, columns C:X) ----
read_header_row <- function(rng) {
  read_excel(pest_data_raw, range = rng, col_names = FALSE, .name_repair = "unique_quiet") |>
    unlist(use.names = FALSE) |>
    as.character()
}
pest_info <- tibble(
  compound = read_header_row("C2:X2"),
  category = read_header_row("C1:X1")
) |>
  fill(category)

pest_names <- pest_info$compound
fung_names <- pest_info$compound[pest_info$category == "Fungicides"]
herb_names <- pest_info$compound[pest_info$category == "Herbicides"]
inse_names <- pest_info$compound[pest_info$category == "Insecticides"]

stopifnot(
  !anyNA(pest_info),
  !anyDuplicated(pest_names),
  length(fung_names) == 14,
  length(herb_names) == 6,
  length(inse_names) == 2
)

# ---- load and tidy ----
# one block of three rows (EtOH, eEtOH, PG) per site and survey date
pest_raw <- read_excel(
  pest_data_raw, skip = 2, col_names = FALSE, .name_repair = "unique_quiet"
) |>
  select(1:(2 + length(pest_names))) |>
  set_names(c("site_raw", "sample_type", pest_names)) |>
  fill(site_raw) |>
  filter(
    str_detect(site_raw, "\\d{2}\\.\\d{2}\\."),
    sample_type %in% names(sample_types)
  ) |>
  mutate(
    Site_ID = str_extract(site_raw, "^[^\n\r]+"),
    survey_date = str_extract(site_raw, "\\d{2}\\.\\d{2}\\."),
    sample_type = unname(sample_types[sample_type]),
    across(all_of(pest_names), as.integer)
  )

stopifnot(
  nrow(pest_raw) == 123,
  all(count(pest_raw, Site_ID, survey_date)$n == 3),
  !anyNA(pest_raw[pest_names]),
  all(unlist(pest_raw[pest_names]) %in% 0:1)
)

# ---- aggregate: union of detections across all survey dates and sample types ----
# n_sampling_dates documents how many survey dates contributed per site
pesticides_all <- pest_raw |>
  group_by(Site_ID) |>
  summarise(
    n_sampling_dates = n_distinct(survey_date),
    across(all_of(pest_names), \(x) as.integer(any(x == 1L))),
    .groups = "drop"
  ) |>
  mutate(
    n_pesticides = as.integer(rowSums(across(all_of(pest_names)))),
    n_fungicides = as.integer(rowSums(across(all_of(fung_names)))),
    n_herbicides = as.integer(rowSums(across(all_of(herb_names)))),
    n_insecticides = as.integer(rowSums(across(all_of(inse_names))))
  ) |>
  select(Site_ID, n_pesticides, n_fungicides, n_herbicides, n_insecticides,
         n_sampling_dates)

# ---- keep breeding sites only and check ----
# MB2 and VB1 were screened but have no confirmed Whinchat breeding record
pesticides <- pesticides_all |>
  filter(Site_ID %in% breeding_ids) |>
  arrange(Site_ID)

stopifnot(
  nrow(pesticides_all) == 22,
  setequal(setdiff(pesticides_all$Site_ID, pesticides$Site_ID), c("MB2", "VB1")),
  nrow(pesticides) == 20,
  setequal(pesticides$Site_ID, breeding_ids),
  !anyDuplicated(pesticides$Site_ID),
  !anyNA(pesticides),
  all(pesticides$n_pesticides ==
        pesticides$n_fungicides + pesticides$n_herbicides + pesticides$n_insecticides)
)

# ---- save ----
saveRDS(pesticides, file.path(output_dir, "pesticides.rds"))
 write_xlsx(pesticides, file.path(output_dir, "pesticides.xlsx"))
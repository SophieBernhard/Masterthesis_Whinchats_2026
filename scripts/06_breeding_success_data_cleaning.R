# 06_breeding_success_data_cleaning.R ----
# Read breeding success data (2024) and aggregate to site level
# 2024 only: all predictors (LiDAR, biomass, metabarcoding, weather,
# vegetation, perches, pesticides) are 2024-specific
# Requires: data_clean/site_ref.rds (produced by script 00)

# ---- packages ----
library("readxl")
library("tidyverse")
library("writexl")
library("here")

# ---- paths ----
breed_data_raw <- here("data_raw/Braunkehlchen_Bruterfolg_2023-24_Auswertung.xlsx")
output_dir <- here("data_clean")

# ---- settings ----
# Raw label for a breeding pair with unresolved outcome (LDK9), the fledgling
# count is set to NA and the site is flagged in fledglings_na
unresolved_label <- "unbekannt, ob erfolgreich"

# ---- name normalisation ----
# identical to the function in script 00
normalise_name <- function(x) {
  x |>
    str_replace_all(c("ü" = "ue", "ö" = "oe", "ä" = "ae", "ß" = "ss")) |>
    str_replace_all("(?<=\\w)-(?=\\w)", " ") |>
    str_replace_all("\\s*:\\s*", ": ") |>
    str_replace_all("\\s+", " ") |>
    str_trim() |>
    str_to_lower()
}

# ---- site reference ----
site_ref <- readRDS(here("data_clean/site_ref.rds"))
breeding_ids <- site_ref |> filter(breeding == "ja") |> pull(Sample_ID)

# ---- raw records 2024 ----
# one row per monitored area, Anzahl = number of pairs or territories
breeding_raw <- read_excel(breed_data_raw, sheet = "Gesamt") |>
  filter(Jahr == 2024) |>
  mutate(
    plot_key = normalise_name(Gebietsname),
    fledglings_raw = str_squish(as.character(`Jungvögel (min.)`)),
    unresolved = fledglings_raw == unresolved_label,
    min_fledglings = as.integer(if_else(unresolved, NA_character_, fledglings_raw)),
    n_pairs = as.integer(Anzahl)
  )

stopifnot(
  nrow(breeding_raw) == 20,
  !anyDuplicated(breeding_raw$plot_key),
  !anyNA(breeding_raw$n_pairs),
  all(breeding_raw$n_pairs >= 1),
  all(is.na(breeding_raw$min_fledglings) == breeding_raw$unresolved),
  all(breeding_raw$min_fledglings >= 0, na.rm = TRUE)
)

# ---- join site IDs ----
breeding_joined <- breeding_raw |>
  left_join(
    site_ref |> filter(breeding == "ja") |> select(plot_key, Sample_ID),
    by = "plot_key",
    relationship = "many-to-one"
  )

stopifnot(!anyNA(breeding_joined$Sample_ID))

# ---- aggregate to site level ----
# total_fledglings is the sum of minimum fledglings across pairs, NA if any
# pair at the site has an unresolved outcome
breeding_success <- breeding_joined |>
  group_by(Site_ID = Sample_ID) |>
  summarise(
    total_fledglings = if (anyNA(min_fledglings)) NA_integer_ else sum(min_fledglings),
    n_pairs_total = sum(n_pairs),
    fledglings_na = anyNA(min_fledglings),
    .groups = "drop"
  ) |>
  arrange(Site_ID)

# ---- check ----
# LDK9 (Driedorf: Mademuehlen Sued 3) is the only site with an unresolved outcome
stopifnot(
  nrow(breeding_success) == 20,
  setequal(breeding_success$Site_ID, breeding_ids),
  !anyDuplicated(breeding_success$Site_ID),
  !anyNA(breeding_success$n_pairs_total),
  setequal(breeding_success$Site_ID[breeding_success$fledglings_na], "LDK9"),
  !anyNA(breeding_success$total_fledglings[!breeding_success$fledglings_na])
)

# ---- save ----
saveRDS(breeding_success, file.path(output_dir, "breeding_success.rds"))
 write_xlsx(breeding_success, file.path(output_dir, "breeding_success.xlsx"))
# 00_site_ref.R ----
# Load, clean, and aggregate site reference data for further use
# in later scripts

# ---- packages ----
library("readxl")
library("tidyverse")
library("here")

# ---- name normalisation ----
normalise_name <- function(x) {
  x |>
    str_replace_all(c("ü" = "ue", "ö"= "oe", "ä" = "ae", "ß" = "ss")) |>
    str_replace_all("(?<=\\w)-(?=\\w)", " ") |>
    str_replace_all("\\s*:\\s*", ": ") |>
    str_replace_all("\\s+", " ") |>
    str_trim() |>
    str_to_lower()
}

# ---- site reference ----
site_ref <- read_excel(here("data_raw/Plot-IDs_Abgleich_alleUntersuchungen.xlsx")) |>
  mutate(
    PlotID_clean = str_remove_all(Flaeche, " "),
    plot_key     = normalise_name(Plot_name)
  ) |>
  select(ID, PlotID_clean, plot_key, Sample_ID, breeding) |>
  distinct()

saveRDS(site_ref, here("data_clean/site_ref.rds"))

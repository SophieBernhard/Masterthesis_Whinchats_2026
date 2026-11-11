# 11_data_aggregation.R ----
# Assemble the two master tables used in all later analysis steps
# master_habitat: n = 40 sites, territory occupation (response breeding_bin)
# master_breeding: n = 19 breeding sites, fledgling output (response total_fledglings, offset n_pairs_total)
# Both tables contain raw predictor values, transformations are applied in script 13
# Requires: data_clean/*.rds produced by scripts 00 to 10

# ---- packages ----
library("tidyverse")
library("writexl")
library("here")

# ---- paths ----
output_dir <- here("data_clean")

# ---- load ----
site_ref <- readRDS(file.path(output_dir, "site_ref.rds"))
veg_diversity <- readRDS(file.path(output_dir, "veg_diversity.rds"))
food_div <- readRDS(file.path(output_dir, "food_diversity_split.rds"))
perches_may <- readRDS(file.path(output_dir, "perches_may.rds"))
perches_breeding <- readRDS(file.path(output_dir, "perches_breeding.rds"))
site_structure <- readRDS(file.path(output_dir, "site_structure.rds"))
biomass <- readRDS(file.path(output_dir, "biomass.rds"))
land_use <- readRDS(file.path(output_dir, "land_use.rds"))
breeding_success <- readRDS(file.path(output_dir, "breeding_success.rds"))
pesticides <- readRDS(file.path(output_dir, "pesticides.rds"))
weather_breeding <- readRDS(file.path(output_dir, "weather_breeding.rds"))

# ---- base tables ----
# ja = confirmed Whinchat breeding site (raw label of the site reference)
base_all <- site_ref |>
  select(Site_ID = Sample_ID, breeding) |>
  distinct() |>
  mutate(
    region = str_extract(Site_ID, "^[A-Z]+"),
    breeding_bin = as.integer(breeding == "ja")
  )

base_breeding <- base_all |> filter(breeding == "ja")

# ---- check (inputs) ----
# every input table holds exactly one row per site and covers all sites it is joined to
has_sites <- function(x, ids) !anyDuplicated(x$Site_ID) && all(ids %in% x$Site_ID)

stopifnot(
  nrow(base_all) == 40,
  !anyDuplicated(base_all$Site_ID),
  nrow(base_breeding) == 20,
  setequal(base_all$region, c("LDK", "MB", "VB")),
  has_sites(veg_diversity, base_all$Site_ID),
  has_sites(food_div, base_all$Site_ID),
  has_sites(perches_may, base_all$Site_ID),
  has_sites(perches_breeding, base_all$Site_ID),
  has_sites(site_structure, base_all$Site_ID),
  has_sites(biomass, base_all$Site_ID),
  has_sites(land_use, base_all$Site_ID),
  has_sites(breeding_success, base_breeding$Site_ID),
  has_sites(pesticides, base_breeding$Site_ID),
  has_sites(weather_breeding, base_breeding$Site_ID)
)

# ---- master_habitat (n = 40) ----
master_habitat <- base_all |>
  select(Site_ID, region, breeding_bin) |>
  left_join(
    food_div |> select(Site_ID, otu_richness_arrival, food_nmds1_arrival,
                       food_nmds2_arrival, food_nmds3_arrival, food_nmds4_arrival,
                       food_beta_mean_arrival),
    by = "Site_ID", relationship = "one-to-one"
  ) |>
  left_join(
    biomass |> select(Site_ID, biomass_mean_arrival),
    by = "Site_ID", relationship = "one-to-one"
  ) |>
  left_join(
    veg_diversity |> select(Site_ID, veg_shannon, veg_richness, veg_evenness),
    by = "Site_ID", relationship = "one-to-one"
  ) |>
  left_join(
    site_structure |> select(Site_ID, viewshed, herb_h_sd, mean_h, max_h,
                             terrain_sd, rumple),
    by = "Site_ID", relationship = "one-to-one"
  ) |>
  left_join(
    perches_may |> select(Site_ID, n_perches),
    by = "Site_ID", relationship = "one-to-one"
  ) |>
  left_join(
    land_use |> select(Site_ID, pct_extensive, pct_intensive, pct_arable,
                       pct_seminat, pct_settlement),
    by = "Site_ID", relationship = "one-to-one"
  ) |>
  arrange(Site_ID)

# ---- master_breeding (n = 19) ----
# LDK9 is dropped because its fledgling count is unresolved (script 06)
master_breeding <- base_breeding |>
  select(Site_ID, region) |>
  left_join(
    breeding_success |> select(Site_ID, total_fledglings, n_pairs_total, fledglings_na),
    by = "Site_ID", relationship = "one-to-one"
  ) |>
  filter(!fledglings_na) |>
  left_join(
    food_div |> select(Site_ID, otu_richness_breeding, food_nmds1_breeding,
                       food_nmds2_breeding, food_nmds3_breeding, food_nmds4_breeding,
                       food_beta_mean_breeding),
    by = "Site_ID", relationship = "one-to-one"
  ) |>
  left_join(
    biomass |> select(Site_ID, biomass_mean_breeding),
    by = "Site_ID", relationship = "one-to-one"
  ) |>
  left_join(
    veg_diversity |> select(Site_ID, veg_shannon, veg_richness, veg_evenness),
    by = "Site_ID", relationship = "one-to-one"
  ) |>
  left_join(
    site_structure |> select(Site_ID, viewshed, herb_h_sd, mean_h, max_h,
                             terrain_sd, rumple, n_perches_25sqm),
    by = "Site_ID", relationship = "one-to-one"
  ) |>
  left_join(
    perches_breeding |> select(Site_ID, n_perches_mid),
    by = "Site_ID", relationship = "one-to-one"
  ) |>
  left_join(
    land_use |> select(Site_ID, pct_extensive, pct_intensive, pct_arable,
                       pct_seminat, pct_settlement),
    by = "Site_ID", relationship = "one-to-one"
  ) |>
  left_join(
    pesticides |> select(Site_ID, n_pesticides, n_fungicides, n_herbicides,
                         n_insecticides),
    by = "Site_ID", relationship = "one-to-one"
  ) |>
  left_join(
    weather_breeding |> select(Site_ID, temp_mean_breeding,
                               max_daily_precip_breeding, n_drydays_breeding,
                               n_heavy_rain_days_breeding),
    by = "Site_ID", relationship = "one-to-one"
  ) |>
  arrange(Site_ID)

# ---- check (result) ----
# LDK17 is the only site without an arrival biomass value (script 08), master_breeding is complete
# variables present in both tables must be identical for the 19 shared sites
shared_vars <- setdiff(intersect(names(master_habitat), names(master_breeding)), "Site_ID")
hab_shared <- master_habitat |>
  filter(Site_ID %in% master_breeding$Site_ID) |>
  select(Site_ID, all_of(shared_vars))

stopifnot(
  nrow(master_habitat) == 40,
  ncol(master_habitat) == 25,
  !anyDuplicated(master_habitat$Site_ID),
  setequal(master_habitat$Site_ID[master_habitat$breeding_bin == 1], base_breeding$Site_ID),
  sum(master_habitat$breeding_bin) == 20,
  setequal(master_habitat$Site_ID[is.na(master_habitat$biomass_mean_arrival)], "LDK17"),
  !anyNA(select(master_habitat, -biomass_mean_arrival)),
  nrow(master_breeding) == 19,
  ncol(master_breeding) == 36,
  !anyDuplicated(master_breeding$Site_ID),
  setequal(master_breeding$Site_ID, setdiff(base_breeding$Site_ID, "LDK9")),
  !anyNA(master_breeding),
  !any(master_breeding$fledglings_na),
  all(master_breeding$n_pairs_total >= 1),
  all(master_breeding$total_fledglings >= 0),
  nrow(hab_shared) == nrow(master_breeding),
  isTRUE(all.equal(as.data.frame(hab_shared),
                   as.data.frame(select(master_breeding, Site_ID, all_of(shared_vars)))))
)

cat("master_habitat: n =", nrow(master_habitat), "| cols =", ncol(master_habitat), "\n")
cat("master_breeding: n =", nrow(master_breeding), "| cols =", ncol(master_breeding), "\n")

# ---- save ----
saveRDS(master_habitat, file.path(output_dir, "master_habitat.rds"))
saveRDS(master_breeding, file.path(output_dir, "master_breeding.rds"))
 write_xlsx(master_habitat, file.path(output_dir, "master_habitat.xlsx"))
 write_xlsx(master_breeding, file.path(output_dir, "master_breeding.xlsx"))
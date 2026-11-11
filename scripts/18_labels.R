# 18_labels.R ----
# Lookup from script codes to figure names
# Save as scripts/labels.R and load with source(here("scripts", "18_labels.R"))

lab <- c(
  # responses
  breeding_bin     = "territory occupation",
  total_fledglings = "fledgling output",
  
  # food PC axes (phase is given by the response, see figure caption)
  food_arrival_PC1  = "food PC1",
  food_arrival_PC2  = "food PC2",
  food_arrival_PC3  = "food PC3",
  food_breeding_PC1 = "food PC1",
  food_breeding_PC2 = "food PC2",
  
  # food variables
  food_nmds1_arrival      = "food NMDS1",
  food_nmds2_arrival      = "food NMDS2",
  food_nmds3_arrival      = "food NMDS3",
  food_nmds4_arrival      = "food NMDS4",
  food_nmds1_breeding     = "food NMDS1",
  food_nmds2_breeding     = "food NMDS2",
  food_nmds3_breeding     = "food NMDS3",
  food_nmds4_breeding     = "food NMDS4",
  biomass_mean_arrival    = "arthropod biomass",
  biomass_mean_breeding   = "arthropod biomass",
  otu_richness_arrival    = "OTU richness",
  otu_richness_breeding   = "OTU richness",
  food_beta_mean_arrival  = "beta diversity",
  food_beta_mean_breeding = "beta diversity",
  
  # vegetation
  veg_shannon  = "veg. Shannon",
  veg_richness = "veg. richness",
  veg_evenness = "veg. evenness",
  
  # habitat structure
  site_structure_arrival_PC1  = "habitat structure PC1",
  site_structure_breeding_PC1 = "habitat structure PC1",
  site_structure_breeding_PC2 = "habitat structure PC2",
  viewshed   = "viewshed",
  herb_h_sd  = "herb SD",
  mean_h     = "mean height",
  max_h      = "max height",
  terrain_sd = "terrain SD",
  rumple     = "rumple",
  
  # perches
  n_perches       = "perch availability",
  n_perches_25sqm = "perches (25 m\u00b2)",
  n_perches_mid   = "perches (mid)",
  
  # land-use
  land_use_arrival_PC1  = "land-use PC1",
  land_use_arrival_PC2  = "land-use PC2",
  land_use_arrival_PC3  = "land-use PC3",
  land_use_breeding_PC1 = "land-use PC1",
  land_use_breeding_PC2 = "land-use PC2",
  pct_extensive = "extensive grassland",
  pct_intensive = "intensive grassland",
  pct_arable    = "arable land",
  pct_seminat   = "semi-natural",
  pct_settlement = "settlement",
  
  # weather
  weather_PC1 = "weather PC1",
  weather_PC2 = "weather PC2",
  temp_mean_breeding        = "mean temperature",
  max_daily_precip_breeding = "max daily precipitation",
  n_drydays_breeding        = "number dry days",
  n_heavy_rain_days_breeding = "number heavy-rain days",
  
  # pesticides
  pesticides_PC1 = "pesticides PC1",
  n_pesticides   = "pesticide",
  n_fungicides   = "fungicide",
  n_herbicides   = "herbicide",
  n_insecticides = "insecticide"
)

# returns figure name, falls back to the code if none is defined
lbl <- function(x) {
  out <- unname(lab[x])
  ifelse(is.na(out), x, out)
}

# DOT node label with line breaks
node_lab <- function(code, width = 16) {
  gsub("\n", "\\\\n", stringr::str_wrap(lbl(code), width))
}
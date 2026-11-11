# 14_PCA_within_variable_groups.R ----
# PCA per predictor group (standardised), retained components: eigenvalue >= 1 (Kaiser)
# Bartlett test of sphericity is reported per group as a diagnostic only, PC signs are arbitrary
# Requires: data_clean/master_habitat_tf.rds and master_breeding_tf.rds (script 13)
# Output: data_clean/pca_scores.rds (read by script 15), output/pca_loadings.xlsx

# ---- packages ----
library("tidyverse")
library("writexl")
library("here")

# ---- paths ----
data_dir <- here("data_clean")
output_dir <- here("output")

# ---- load ----
hab <- readRDS(file.path(data_dir, "master_habitat_tf.rds"))
bre <- readRDS(file.path(data_dir, "master_breeding_tf.rds"))
stopifnot(!anyDuplicated(hab$Site_ID), !anyDuplicated(bre$Site_ID))

# ---- define groups ----
groups <- list(
  
  # ---- habitat (arrival) ----
  list(src = hab, group = "food_arrival",
       vars = c("otu_richness_arrival", "food_nmds1_arrival", "food_nmds2_arrival",
                "food_nmds3_arrival", "food_nmds4_arrival",
                "food_beta_mean_arrival", "biomass_mean_arrival")),
  
  list(src = hab, group = "site_structure_arrival",
       vars = c("viewshed", "herb_h_sd", "mean_h_log", "terrain_sd_log", "rumple_log")),
  
  list(src = hab, group = "land_use_arrival",
       vars = c("pct_extensive", "pct_intensive", "pct_arable",
                "pct_seminat_log", "pct_settlement_log")),
  
  # ---- breeding ----
  list(src = bre, group = "food_breeding",
       vars = c("otu_richness_breeding", "food_nmds1_breeding", "food_nmds2_breeding",
                "food_nmds3_breeding", "food_nmds4_breeding",
                "food_beta_mean_breeding", "biomass_mean_breeding_log")),
  
  list(src = bre, group = "site_structure_breeding",
       vars = c("viewshed", "herb_h_sd", "mean_h_log", "terrain_sd", "rumple")),
  
  list(src = bre, group = "land_use_breeding",
       vars = c("pct_extensive", "pct_intensive", "pct_arable",
                "pct_seminat_log", "pct_settlement_log")),
  
  list(src = bre, group = "pesticides",
       vars = c("n_fungicides", "n_herbicides", "n_insecticides")),
  
  list(src = bre, group = "weather",
       vars = c("temp_mean_breeding", "max_daily_precip_breeding",
                "n_drydays_breeding", "n_heavy_rain_days_breeding"))
)

# ---- Bartlett test ----
bartlett <- function(mat) {
  n <- nrow(mat); p <- ncol(mat)
  chi <- -(n - 1 - (2 * p + 5) / 6) * log(det(cor(mat)))
  df <- p * (p - 1) / 2
  tibble(n = n, p = p, n_per_var = round(n / p, 1),
         bartlett_chi2 = round(chi, 2), bartlett_df = df,
         bartlett_p = signif(pchisq(chi, df, lower.tail = FALSE), 3))
}

# ---- run PCA ----
run_pca <- function(g) {
  tmp <- g$src |> select(Site_ID, all_of(g$vars)) |> na.omit()
  mat <- tmp |> select(all_of(g$vars)) |> as.matrix()
  pca <- prcomp(mat, scale. = TRUE, center = TRUE)
  eval <- pca$sdev^2
  npcs <- sum(eval >= 1)
  
  loads <- as.data.frame(pca$rotation[, seq_len(npcs), drop = FALSE]) |>
    rownames_to_column("variable") |>
    mutate(group = g$group, .before = 1)
  
  scores <- as.data.frame(pca$x[, seq_len(npcs), drop = FALSE]) |>
    set_names(paste0(g$group, "_PC", seq_len(npcs))) |>
    mutate(Site_ID = tmp$Site_ID, .before = 1)
  
  var_tab <- tibble(group = g$group, component = paste0("PC", seq_along(eval)),
                    eigenvalue = round(eval, 3),
                    var_pct = round(eval / sum(eval) * 100, 1),
                    cum_pct = round(cumsum(eval) / sum(eval) * 100, 1),
                    retained = seq_along(eval) <= npcs)
  
  list(group = g$group, loads = loads, scores = scores, var_tab = var_tab,
       diag = bind_cols(tibble(group = g$group), bartlett(mat)))
}

pca_list <- map(groups, run_pca)

# ---- save ----
saveRDS(set_names(map(pca_list, "scores"), map_chr(pca_list, "group")),
        file.path(data_dir, "pca_scores.rds"))

write_xlsx(set_names(map(pca_list, "scores"), map_chr(pca_list, "group")),
        file.path(output_dir, "pca_scores.xlsx"))
write_xlsx(list(loadings = bind_rows(map(pca_list, "loads")),
                variance = bind_rows(map(pca_list, "var_tab")),
                diagnostics = bind_rows(map(pca_list, "diag"))),
           file.path(output_dir, "pca_loadings.xlsx"))
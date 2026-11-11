# 15_cross_group_screening.R ----
# Cross-group Spearman screening, VIF and a-priori reduction of the PCA-based predictors
# Requires: master_habitat_tf.rds, master_breeding_tf.rds (script 13), pca_scores.rds (script 14)
# Output: pred_hab.rds, pred_bre.rds (scripts 16 and 17), corr_cross_groups.xlsx, vif.xlsx

# ---- packages ----
library("tidyverse")
library("writexl")
library("here")

# ---- settings ----
data_dir <- here("data_clean")
output_dir <- here("output")
rho_threshold <- 0.7
vif_max <- 5
excl_hab <- "LDK17" # only site without food arrival scores, excluded from all habitat analyses

# a-priori exclusions (step 2 low explained variance, step 3 ecological redundancy)
drop_var <- list(hab = c("food_arrival_PC3", "land_use_arrival_PC3"),
                 bre = c("site_structure_breeding_PC2", "land_use_breeding_PC2",
                         "land_use_breeding_PC3", "food_breeding_PC2"))
drop_eco <- list(hab = character(0),
                 bre = c("site_structure_breeding_PC1", "land_use_breeding_PC1"))

# ---- predictor tables ----
hab <- readRDS(file.path(data_dir, "master_habitat_tf.rds"))
bre <- readRDS(file.path(data_dir, "master_breeding_tf.rds"))
scores <- readRDS(file.path(data_dir, "pca_scores.rds"))

pred_hab <- hab |>
  select(Site_ID, region, breeding_bin, veg_shannon, n_perches) |>
  left_join(reduce(scores[c("food_arrival", "site_structure_arrival", "land_use_arrival")],
                   full_join, by = "Site_ID"), by = "Site_ID")

pred_bre <- bre |>
  select(Site_ID, region, total_fledglings, n_pairs_total, veg_shannon, n_perches_25sqm) |>
  left_join(reduce(scores[c("food_breeding", "site_structure_breeding", "land_use_breeding",
                            "pesticides", "weather")], full_join, by = "Site_ID"), by = "Site_ID")

preds_hab <- setdiff(names(pred_hab), c("Site_ID", "region", "breeding_bin"))
preds_bre <- setdiff(names(pred_bre), c("Site_ID", "region", "total_fledglings", "n_pairs_total"))

# the joins must not change the number of sites, and LDK17 must be the only incomplete site
stopifnot(nrow(pred_hab) == nrow(hab), nrow(pred_bre) == nrow(bre), !anyNA(pred_bre), !anyNA(pred_hab$breeding_bin),
          identical(pred_hab$Site_ID[!complete.cases(pred_hab[preds_hab])], excl_hab),
          all(c(unlist(drop_var), unlist(drop_eco)) %in% c(preds_hab, preds_bre)))

# LDK17 is removed before screening, so pairs, VIF and models in scripts 16 and 17 share n = 39
pred_hab <- filter(pred_hab, Site_ID != excl_hab)

# ---- cross-group Spearman ----
# BH-FDR is applied across all pairs within a dataset
cross_corr <- function(df, preds) {
  grp <- set_names(sub("_PC[0-9]+$", "", preds), preds)
  pairs <- Filter(\(p) grp[[p[1]]] != grp[[p[2]]], combn(preds, 2, simplify = FALSE))
  map(pairs, \(p) {
    tst <- cor.test(df[[p[1]]], df[[p[2]]], method = "spearman", exact = FALSE)
    tibble(var1 = p[1], var2 = p[2], group1 = grp[[p[1]]], group2 = grp[[p[2]]],
           n = nrow(df), rho = unname(tst$estimate), p_raw = tst$p.value)
  }) |>
    list_rbind() |>
    mutate(p_fdr = p.adjust(p_raw, "BH")) |>
    arrange(desc(abs(rho))) |>
    mutate(flag = abs(rho) >= rho_threshold,
           rho = round(rho, 3), p_raw = round(p_raw, 4), p_fdr = round(p_fdr, 4))
}

corr_hab <- cross_corr(pred_hab, preds_hab)
corr_bre <- cross_corr(pred_bre, preds_bre)
stopifnot(nrow(corr_hab) == 30, nrow(corr_bre) == 60, !any(corr_hab$flag), !any(corr_bre$flag))

# ---- VIF ----
vif_cormat <- function(df, preds) {
  tibble(predictor = preds, VIF = diag(solve(cor(df[preds])))) |>
    arrange(desc(VIF)) |>
    mutate(flag = VIF > vif_max, VIF = round(VIF, 2))
}

vif_hab <- vif_cormat(pred_hab, preds_hab)
vif_bre <- vif_cormat(pred_bre, preds_bre)

# ---- a-priori reduction ----
# step 1 VIF filter, step 2 low explained variance, step 3 ecological redundancy
keep_vif <- list(hab = vif_hab$predictor[!vif_hab$flag], bre = vif_bre$predictor[!vif_bre$flag])
keep_fin <- map2(map2(keep_vif, drop_var, setdiff), drop_eco, setdiff)

# ---- save ----
saveRDS(pred_hab, file.path(output_dir, "pred_hab.rds"))
saveRDS(pred_bre, file.path(output_dir, "pred_bre.rds"))
write_xlsx(list(habitat = corr_hab, breeding = corr_bre),
           file.path(output_dir, "corr_cross_groups.xlsx"))
write_xlsx(list(habitat = vif_hab, breeding = vif_bre),
           file.path(output_dir, "vif.xlsx"))
write_xlsx(pred_hab, file.path(output_dir, "pred_hab.xlsx"))
write_xlsx(pred_bre, file.path(output_dir, "pred_bre.xlsx"))

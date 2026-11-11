# 12_within_group_correlations.R ----
# Pairwise Spearman correlations within predictor groups (collinearity screening)
# Justifies the variable exclusions (veg_richness, veg_evenness, max_h, n_pesticides)
# and the group PCAs in script 14
# Requires: data_clean/master_habitat.rds and master_breeding.rds (script 11),
#           scripts/18_labels.R (lbl())
# Spearman is rank-based, so results are identical before and after the monotonic
# transformations of script 13. Ties are expected (pct_* zeros, binary n_insecticides,
# few-valued n_perches_mid and n_heavy_rain_days), hence exact = FALSE (normal approximation)
# BH correction is applied within each panel (= correlation matrix)

# ---- packages ----
library("tidyverse")
library("here")
library("ggcorrplot")
library("writexl")
library("patchwork")
library("ragg")

# ---- paths ----
data_dir <- here("data_clean")
output_dir <- here("output")
source(here("scripts", "18_labels.R"))

# ---- settings ----
rho_threshold <- 0.7   # |rho| >= threshold flags a collinear pair
sig_level <- 0.05      # BH-adjusted p used to show cells in the plots
n_col <- 4             # panel columns in the combined figure

# ---- load ----
df_h <- readRDS(file.path(data_dir, "master_habitat.rds"))
df_b <- readRDS(file.path(data_dir, "master_breeding.rds"))

# ---- variable sets used in both datasets ----
v_veg <- c("veg_shannon", "veg_richness", "veg_evenness")
v_site <- c("viewshed", "herb_h_sd", "mean_h", "max_h", "terrain_sd", "rumple")
v_land <- c("pct_extensive", "pct_intensive", "pct_arable", "pct_seminat", "pct_settlement")

# ---- panels = tested groups per dataset (order = plot order) ----
# row 1 = territory occupation (df_h, n = 40), rows 2-3 = fledgling output (df_b, n = 19)
# food order as in scripts 13 and 14
# n_pesticides is the sum of the three pesticide counts (collinear by construction)
panels <- list(
  food_arrival = list(df = df_h, title = "food availability territory occupation",
                      vars = c("otu_richness_arrival", "food_nmds1_arrival",
                               "food_nmds2_arrival", "food_nmds3_arrival",
                               "food_nmds4_arrival", "food_beta_mean_arrival",
                               "biomass_mean_arrival")),
  veg_arrival = list(df = df_h, vars = v_veg, title = "vegetation territory occupation"),
  site_structure_arrival = list(df = df_h, vars = v_site,
                                title = "habitat structure territory occupation"),
  land_use_arrival = list(df = df_h, vars = v_land, title = "land use territory occupation"),
  food_breeding = list(df = df_b, title = "food availability fledgling output",
                       vars = c("otu_richness_breeding", "food_nmds1_breeding",
                                "food_nmds2_breeding", "food_nmds3_breeding",
                                "food_nmds4_breeding", "food_beta_mean_breeding",
                                "biomass_mean_breeding")),
  veg_breeding = list(df = df_b, vars = v_veg, title = "vegetation fledgling output"),
  site_structure_breeding = list(df = df_b, vars = v_site,
                                 title = "habitat structure fledgling output"),
  land_use_breeding = list(df = df_b, vars = v_land, title = "land use fledgling output"),
  perches_breeding = list(df = df_b, vars = c("n_perches_25sqm", "n_perches_mid"),
                          title = "perch availability fledgling output"),
  weather = list(df = df_b, title = "weather fledgling output",
                 vars = c("temp_mean_breeding", "max_daily_precip_breeding",
                          "n_drydays_breeding", "n_heavy_rain_days_breeding")),
  pesticides = list(df = df_b, title = "pesticides fledgling output",
                    vars = c("n_pesticides", "n_fungicides", "n_herbicides", "n_insecticides"))
)

# ---- check ----
stopifnot(
  nrow(df_h) == 40, nrow(df_b) == 19,
  !anyNA(df_h$Site_ID), !anyNA(df_b$Site_ID),
  !anyDuplicated(df_h$Site_ID), !anyDuplicated(df_b$Site_ID),
  length(panels) == 11
)
for (nm in names(panels)) {
  p <- panels[[nm]]
  stopifnot(
    all(p$vars %in% names(p$df)),                       # no silent dropping of variables
    all(map_lgl(p$df[p$vars], is.numeric)),
    all(map_int(p$df[p$vars], \(x) n_distinct(x, na.rm = TRUE)) > 1),  # no constant variable
    !anyNA(lbl(p$vars))                                 # label exists for every variable
  )
}
# only expected NA: biomass_mean_arrival at LDK17 (n = 39 in biomass pairs)
na_vars <- unique(unlist(map(panels, \(p) p$vars[map_lgl(p$df[p$vars], anyNA)])))
stopifnot(identical(na_vars, "biomass_mean_arrival"))

# ---- pairwise Spearman + BH correction within each panel ----
corr_pairs <- function(df, vars, label) {
  combn(vars, 2, simplify = FALSE) |>
    map(\(p) {
      cc <- complete.cases(df[p])
      stopifnot(sum(cc) >= 4)
      tst <- cor.test(df[[p[1]]][cc], df[[p[2]]][cc], method = "spearman", exact = FALSE)
      tibble(group = label, var1 = p[1], var2 = p[2], n = sum(cc),
             rho = unname(tst$estimate), p_raw = tst$p.value)
    }) |>
    list_rbind() |>
    mutate(
      p_fdr = p.adjust(p_raw, "BH"),
      flag_collinear = abs(rho) >= rho_threshold,   # flag from unrounded rho
      sig = as.character(symnum(p_fdr, cutpoints = c(0, .001, .01, .05, .1, 1),
                                symbols = c("***", "**", "*", ".", "ns"), legend = FALSE))
    )
}

results <- imap(panels, \(p, nm) corr_pairs(p$df, p$vars, nm)) |>
  list_rbind() |>
  mutate(rho = round(rho, 3))   # round only for reporting, p-values stay unrounded

collinear <- filter(results, flag_collinear)
stopifnot(nrow(results) == 111, nrow(collinear) == 8)

message("\nCollinear pairs (|rho| >= ", rho_threshold, "):")
print(select(collinear, group, var1, var2, n, rho, p_fdr, sig), n = 50)
print(results, n = 150)

# ---- correlation matrix plots built from results (BH-adjusted p) ----
corrmat_plot <- function(p, nm) {
  res <- filter(results, group == nm)
  vars <- p$vars
  rho <- diag(length(vars))
  pmat <- matrix(1, length(vars), length(vars))
  diag(pmat) <- 0
  idx <- cbind(match(res$var1, vars), match(res$var2, vars))
  rho[idx] <- rho[idx[, 2:1]] <- res$rho
  pmat[idx] <- pmat[idx[, 2:1]] <- res$p_fdr
  dimnames(rho) <- dimnames(pmat) <- list(lbl(vars), lbl(vars))
  ggcorrplot(
    rho, p.mat = pmat, type = "lower", lab = TRUE, lab_size = 4.2,
    insig = "blank", sig.level = sig_level,
    colors = c("#D55E00", "white", "#009E73"),
    ggtheme = theme_minimal(base_size = 10, base_family = "Calibri")
  ) +
    labs(caption = p$title) +
    theme(
      plot.caption = element_text(face = "bold", size = 12, family = "Calibri", hjust = 0.5),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 12, family = "Calibri"),
      axis.text.y = element_text(size = 12, family = "Calibri"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank()
    )
}

plots <- imap(panels, corrmat_plot)
stopifnot(length(plots) == 11)

# ---- save single panels ----
single_dir <- file.path(output_dir, "corr_within_groups")
dir.create(single_dir, showWarnings = FALSE)

iwalk(plots, \(pl, nm) {
  k <- length(panels[[nm]]$vars)
  ggsave(file.path(single_dir, paste0("corr_", nm, ".png")),
         plot = pl, width = 2.5 + 0.55 * k, height = 2.5 + 0.55 * k,
         device = ragg::agg_png, dpi = 300)
})

# ---- save ----
 write_xlsx(results, file.path(output_dir, "corr_within_groups.xlsx"))
# 16_model_selection_dredge.R ----
# Dredge and model averaging for territory occupation (binomial GLM) and fledgling output (Poisson GLM)
# Requires: output/pred_hab.rds, pred_bre.rds (script 15), data_raw coordinates (mastertable.xlsx)
# Output: dredge_hab.xlsx, dredge_bre.xlsx, diagnostics_16.xlsx, diag_residuals.png, diag_roc_hab.png, fig_marginal_effects_hab.png

# ---- packages ----
library("tidyverse")
library("here")
library("writexl")
library("MuMIn")
library("ggeffects")
library("patchwork")
library("pROC")
library("logistf")
library("ResourceSelection")
library("ragg")

# ---- paths ----
output_dir <- here("output")

# ---- settings ----
delta_cut <- 4 # models within this delta AICc form the competitive set
hl_groups <- 10 # groups for the Hosmer-Lemeshow test
z_crit <- qnorm(0.975) # critical value for Wald-type 95% intervals
exp_n <- c(hab = 39L, bre = 19L) # expected number of sites per model
exp_occupied <- 19L # expected number of occupied habitat sites
exp_models <- c(hab = 128L, bre = 64L) # 2^7 and 2^6 candidate models
excl_sens <- "MB7" # site with 10 fledglings from 3 pairs, excluded in the sensitivity check (plus highest-leverage site)
hab_preds <- c("food_arrival_PC1", "food_arrival_PC2", "site_structure_arrival_PC1",
               "land_use_arrival_PC1", "land_use_arrival_PC2", "veg_shannon", "n_perches")
bre_preds <- c("food_breeding_PC1", "n_perches_25sqm", "pesticides_PC1",
               "weather_PC1", "weather_PC2", "veg_shannon")
fig_preds <- c("n_perches", "food_arrival_PC2", "site_structure_arrival_PC1") # possible figure predictors, in panel order
col_hab <- "#009E73"
col_smooth <- "#D55E00"

# ---- plot theme ----
theme_thesis <- function(size = 10) {
  theme_bw(base_size = size, base_family = "Calibri") +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = size, family = "Calibri"),
          plot.title = element_blank(),
          axis.title = element_text(size = size, family = "Calibri"),
          axis.title.x = element_text(margin = margin(t = 8)),
          axis.title.y = element_text(margin = margin(r = 8)),
          axis.text = element_text(size = size, family = "Calibri"))
}

# ---- load ----
pred_hab <- readRDS(file.path(output_dir, "pred_hab.rds"))
pred_bre <- readRDS(file.path(output_dir, "pred_bre.rds"))

source(here("scripts", "19_helper_coordinates.R"))
source(here("scripts", "18_labels.R"))

options(na.action = na.fail)

# ---- input checks ----
stopifnot(nrow(pred_hab) == exp_n[["hab"]], nrow(pred_bre) == exp_n[["bre"]],
          sum(pred_hab$breeding_bin) == exp_occupied,
          !anyNA(pred_hab[hab_preds]), !anyNA(pred_bre[bre_preds]),
          !anyDuplicated(coords_df$Site_ID),
          all(c(pred_hab$Site_ID, pred_bre$Site_ID) %in% coords_df$Site_ID),
          !anyNA(coords_df[coords_df$Site_ID %in% c(pred_hab$Site_ID, pred_bre$Site_ID), ]),
          excl_sens %in% pred_bre$Site_ID)

# ---- helper functions ----
coef_table <- function(avg) {
  summary(avg)$coefmat.full |>
    as.data.frame() |>
    setNames(c("estimate", "se", "adj_se", "z", "p")) |>
    rownames_to_column("predictor") |>
    filter(predictor != "(Intercept)") |>
    mutate(ci_low = estimate - z_crit * adj_se,
           ci_high = estimate + z_crit * adj_se) |>
    arrange(desc(abs(estimate)))
}

imp_table <- function(avg) {
  tibble(predictor = names(sw(avg)), importance = as.numeric(sw(avg))) |>
    filter(!grepl("^offset", predictor))
}

infl_summary <- function(mod) {
  tibble(max_cooks_d = max(cooks.distance(mod)),
         max_leverage = max(hatvalues(mod)),
         leverage_threshold = 2 * length(coef(mod)) / nobs(mod))
}

# ---- habitat model (binomial GLM, breeding_bin) ----
dat_hab_full <- pred_hab |>
  select(Site_ID, region, breeding_bin, all_of(hab_preds)) |>
  mutate(across(all_of(hab_preds), ~as.numeric(scale(.))))

dat_hab <- select(dat_hab_full, -Site_ID, -region)

glm_hab <- glm(breeding_bin ~ ., data = dat_hab, family = binomial)
dr_hab <- dredge(glm_hab)
top_hab <- subset(dr_hab, delta < delta_cut)
avg_hab <- model.avg(top_hab)
summary(avg_hab)

coef_hab <- coef_table(avg_hab)
imp_hab <- imp_table(avg_hab)

# ---- breeding model (Poisson with offset, total_fledglings) ----
dat_bre_full <- pred_bre |>
  select(Site_ID, region, total_fledglings, n_pairs_total, all_of(bre_preds)) |>
  mutate(across(all_of(bre_preds), ~as.numeric(scale(.))))

dat_bre <- dat_bre_full |>
  select(-Site_ID, -region) |>
  mutate(log_n_pairs = log(n_pairs_total))

glm_bre <- glm(total_fledglings ~ . - n_pairs_total - log_n_pairs + offset(log_n_pairs),
               data = dat_bre, family = poisson)
c_hat <- sum(residuals(glm_bre, type = "pearson")^2) / df.residual(glm_bre)

dr_bre <- dredge(glm_bre, fixed = "offset(log_n_pairs)")
top_bre <- subset(dr_bre, delta < delta_cut)
avg_bre <- model.avg(top_bre)
summary(avg_bre)

coef_bre <- coef_table(avg_bre)
imp_bre <- imp_table(avg_bre)

# ---- best models ----
stopifnot(nrow(dr_hab) == exp_models[["hab"]], nrow(dr_bre) == exp_models[["bre"]])
best_hab <- get.models(dr_hab, subset = 1)[[1]]
best_bre <- get.models(dr_bre, subset = 1)[[1]]

# ---- goodness of fit and influence ----
hl_hab <- hoslem.test(dat_hab$breeding_bin, fitted(glm_hab), g = hl_groups)
print(hl_hab)
cat("c-hat (global fledgling model):", round(c_hat, 3), "\n")

infl_hab <- infl_summary(best_hab)
infl_bre <- infl_summary(best_bre)
print(infl_hab)
print(infl_bre)

# ---- Firth penalized regression, sensitivity check (habitat model only) ----
firth_hab <- logistf(formula(best_hab), data = dat_hab)

firth_compare <- tibble(
  predictor = names(coef(best_hab)),
  glm_coef = round(coef(best_hab), 4),
  firth_coef = round(coef(firth_hab), 4),
  abs_diff = round(abs(coef(best_hab) - coef(firth_hab)), 4)
)
print(firth_compare)

# ---- sensitivity checks without influential sites (breeding model) ----
# MB7 (high output) and the site with the highest leverage are removed singly and together
lev_site <- dat_bre_full$Site_ID[which.max(hatvalues(best_bre))]
excl_all <- unique(c(excl_sens, lev_site))
sens_sets <- c(setNames(as.list(excl_all), excl_all), list(both = excl_all))
coef_row <- function(mod, label) {
  tibble(excluded = label, predictor = names(coef(mod)), coef = coef(mod),
         p = summary(mod)$coefficients[, 4])
}
sens_bre <- bind_rows(
  coef_row(best_bre, "none"),
  list_rbind(imap(sens_sets, \(ex, nm) coef_row(glm(formula(best_bre), family = poisson,
                                                    data = dat_bre[!dat_bre_full$Site_ID %in% ex, ]), nm)))) |>
  mutate(across(where(is.numeric), ~round(., 4)))
cat("Highest-leverage site:", lev_site, "\n")
print(sens_bre)

# ---- residual plots ----
make_resid_plot <- function(mod_label, fitted_vals, resid_vals) {
  tibble(panel = mod_label, fitted = fitted_vals, residual = resid_vals) |>
    ggplot(aes(fitted, residual)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = col_smooth,
                linewidth = 0.8, span = 1) +
    facet_wrap(~panel) +
    xlab("fitted values") + ylab("deviance residuals") +
    theme_thesis(9)
}

p_resid <- make_resid_plot(sprintf("territory occupation (n = %d)", nrow(dat_hab)),
                           fitted(best_hab), residuals(best_hab, type = "deviance")) +
  make_resid_plot(sprintf("fledgling output (n = %d)", nrow(dat_bre)),
                  fitted(best_bre), residuals(best_bre, type = "deviance"))

# ---- ROC / AUC ----
roc_hab <- roc(dat_hab$breeding_bin, fitted(best_hab), quiet = TRUE)
auc_hab <- as.numeric(auc(roc_hab))
ci_hab <- as.numeric(ci.auc(roc_hab, method = "delong"))[c(1, 3)] # lower, upper
cat("Training AUC:", round(auc_hab, 3), "| DeLong 95% CI:", round(ci_hab, 3), "\n")

p_roc <- ggroc(roc_hab, colour = col_hab, linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", colour = "grey55") +
  xlab("specificity") + ylab("sensitivity") +
  theme_thesis(9)

# ---- LOO-CV AUC ----
# fixed formula (best model of the full data) and nested (all subsets refitted in every fold)
loo_fixed <- vapply(seq_len(nrow(dat_hab)), function(i) {
  m_i <- glm(formula(best_hab), data = dat_hab[-i, ], family = binomial)
  as.numeric(predict(m_i, newdata = dat_hab[i, ], type = "response"))
}, numeric(1))

subsets_hab <- unlist(lapply(0:length(hab_preds), \(k) combn(hab_preds, k, simplify = FALSE)),
                      recursive = FALSE)

loo_nested <- vapply(seq_len(nrow(dat_hab)), function(i) {
  dat_tr <- dat_hab[-i, ]
  fits <- lapply(subsets_hab, \(v) glm(reformulate(if (length(v)) v else "1", "breeding_bin"),
                                       data = dat_tr, family = binomial))
  m_i <- fits[[which.min(vapply(fits, AICc, numeric(1)))]]
  as.numeric(predict(m_i, newdata = dat_hab[i, ], type = "response"))
}, numeric(1))

auc_loo_fixed <- as.numeric(auc(roc(dat_hab$breeding_bin, loo_fixed, quiet = TRUE)))
auc_loo_nested <- as.numeric(auc(roc(dat_hab$breeding_bin, loo_nested, quiet = TRUE)))
cat("LOO-CV AUC, fixed formula:", round(auc_loo_fixed, 3), "\n")
cat("LOO-CV AUC, nested selection:", round(auc_loo_nested, 3), "\n")

# ---- spatial autocorrelation (Moran's I) ----
Moran_hab <- moran_resid(best_hab, dat_hab_full$Site_ID)
Moran_bre <- moran_resid(best_bre, dat_bre_full$Site_ID)
print(Moran_hab)
print(Moran_bre)

# ---- region test (Kruskal-Wallis) ----
kw_hab <- kruskal.test(residuals(best_hab, type = "deviance") ~ dat_hab_full$region)
kw_bre <- kruskal.test(residuals(best_bre, type = "deviance") ~ dat_bre_full$region)
print(kw_hab)
print(kw_bre)

# ---- marginal effects figure, territory occupation (best model) ----
fig_vars <- intersect(fig_preds, names(coef(best_hab)))
stopifnot(setequal(fig_vars, setdiff(names(coef(best_hab)), "(Intercept)")))

p_fig <- wrap_plots(imap(set_names(fig_vars), \(v, nm) {
  plot(ggpredict(best_hab, terms = paste0(v, " [all]")), show_data = TRUE, colors = col_hab) +
    coord_cartesian(ylim = c(0, 1)) +
    xlab(paste0(lbl(v), " (z-score)")) + ylab("predicted probability") +
    theme_thesis(9)
}), ncol = length(fig_vars))

# ---- key results table ----
key_results <- bind_rows(
  bind_cols(
    tibble(model = "territory occupation", n_sites = nrow(dat_hab), n_competitive = nrow(top_hab),
           hl_chisq = unname(hl_hab$statistic), hl_p = hl_hab$p.value,
           auc_training = auc_hab, auc_ci_low = ci_hab[1], auc_ci_high = ci_hab[2],
           auc_loo_fixed, auc_loo_nested),
    infl_hab,
    tibble(moran_I = Moran_hab$observed, moran_p = Moran_hab$p.value, kw_p = kw_hab$p.value)),
  bind_cols(
    tibble(model = "fledgling output", n_sites = nrow(dat_bre), n_competitive = nrow(top_bre),
           c_hat = c_hat),
    infl_bre,
    tibble(moran_I = Moran_bre$observed, moran_p = Moran_bre$p.value, kw_p = kw_bre$p.value)))
print(key_results)

# ---- save ----
ggsave(file.path(output_dir, "diag_residuals.png"), p_resid,
       width = 6.5, height = 3, device = ragg::agg_png, dpi = 300)
ggsave(file.path(output_dir, "diag_roc_hab.png"), p_roc,
       width = 4, height = 4, device = ragg::agg_png, dpi = 300)
ggsave(file.path(output_dir, "fig_marginal_effects_hab.png"), p_fig,
       width = 5, height = 3, device = ragg::agg_png, dpi = 300)

write_xlsx(list(all_models = as.data.frame(dr_hab),
 top_delta4 = as.data.frame(top_hab),
 coef_avg = coef_hab,
 importance = imp_hab),
 file.path(output_dir, "dredge_hab.xlsx"))

write_xlsx(list(all_models = as.data.frame(dr_bre),
 top_delta4 = as.data.frame(top_bre),
 coef_avg = coef_bre,
 importance = imp_bre),
 file.path(output_dir, "dredge_bre.xlsx"))

write_xlsx(list(key_results = key_results, firth = firth_compare,
 sensitivity_bre = sens_bre),
 file.path(output_dir, "diagnostics_16.xlsx"))
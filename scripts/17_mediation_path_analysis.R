# 17_mediation_path_analysis.R ----
# Requires: output/pred_hab.rds, pred_bre.rds (script 15), scripts/19_helper_coordinates.R
# Output: sem_path_hab.png, sem_path_bre.png, qqplots_sem_resid.png, sem_results.xlsx

library("tidyverse")
library("here")
library("piecewiseSEM")
library("lmtest")
library("writexl")
library("DiagrammeR")
library("DiagrammeRsvg")
library("rsvg")
library("ragg")

source(here("scripts", "19_helper_coordinates.R"))
source(here("scripts", "18_labels.R"))

output_dir <- here("output")
alpha <- 0.05
n_boot <- 999L
col_resp <- "#009E73"
col_food <- "#D55E00"
fill_food <- "#FFE8D6"

theme_thesis <- function(size = 10) {
  theme_bw(base_size = size, base_family = "Calibri") +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = size), plot.title = element_blank(),
          axis.text = element_text(size = size),
          axis.title.x = element_text(margin = margin(t = 8)),
          axis.title.y = element_text(margin = margin(r = 8)))
}

# ---- data ----
hab_vars <- c("Site_ID", "breeding_bin", "n_perches", "site_structure_arrival_PC1",
              "food_arrival_PC1", "food_arrival_PC2", "land_use_arrival_PC1",
              "land_use_arrival_PC2", "veg_shannon")
bre_vars <- c("Site_ID", "total_fledglings", "n_pairs_total", "food_breeding_PC1",
              "pesticides_PC1", "n_perches_25sqm", "weather_PC1", "weather_PC2", "veg_shannon")

pred_hab <- readRDS(file.path(output_dir, "pred_hab.rds")) |> select(all_of(hab_vars))
pred_bre <- readRDS(file.path(output_dir, "pred_bre.rds")) |>
  select(all_of(bre_vars)) |>
  mutate(log_n_pairs = log(n_pairs_total))

stopifnot(nrow(pred_hab) == 39, nrow(pred_bre) == 19, !anyNA(pred_hab), !anyNA(pred_bre),
          all(c(pred_hab$Site_ID, pred_bre$Site_ID) %in% coords_df$Site_ID))

# ---- models ----
m_bre_resp <- glm(total_fledglings ~ food_breeding_PC1 + pesticides_PC1 + n_perches_25sqm +
                    weather_PC1 + weather_PC2 + veg_shannon + offset(log_n_pairs),
                  family = poisson, data = pred_bre)
m_food_bre1 <- lm(food_breeding_PC1 ~ weather_PC1 + weather_PC2 + pesticides_PC1 + veg_shannon,
                  data = pred_bre)
m_hab_resp <- glm(breeding_bin ~ n_perches + site_structure_arrival_PC1 + food_arrival_PC1 +
                    food_arrival_PC2 + land_use_arrival_PC1 + land_use_arrival_PC2 + veg_shannon,
                  family = binomial, data = pred_hab)
m_food_hab1 <- lm(food_arrival_PC1 ~ site_structure_arrival_PC1 + land_use_arrival_PC1 +
                    land_use_arrival_PC2 + veg_shannon, data = pred_hab)
m_food_hab2 <- lm(food_arrival_PC2 ~ site_structure_arrival_PC1 + land_use_arrival_PC1 +
                    land_use_arrival_PC2 + veg_shannon, data = pred_hab)

sem_bre <- psem(m_bre_resp, m_food_bre1, data = pred_bre)
sem_hab <- psem(m_hab_resp, m_food_hab1, m_food_hab2, data = pred_hab)

name_coefs <- function(x) {
  x <- as.data.frame(x)
  names(x)[names(x) == ""] <- "signif"
  as_tibble(x)
}
coef_bre <- name_coefs(coefs(sem_bre))
coef_hab <- name_coefs(coefs(sem_hab))

indirect_effects <- function(cf, resp, meds) {
  a <- cf |>
    filter(Response %in% meds) |>
    select(predictor = Predictor, mediator = Response, a_std = Std.Estimate)
  b <- cf |>
    filter(Response == resp, Predictor %in% meds) |>
    select(mediator = Predictor, b_std = Std.Estimate)
  inner_join(a, b, by = "mediator") |>
    mutate(response = resp, indirect_std = a_std * b_std) |>
    select(predictor, mediator, response, a_std, b_std, indirect_std)
}
indirect_bre <- indirect_effects(coef_bre, "total_fledglings", "food_breeding_PC1")
indirect_hab <- indirect_effects(coef_hab, "breeding_bin", c("food_arrival_PC1", "food_arrival_PC2"))

r_squared <- bind_rows(mutate(rsquared(sem_bre), sem = "breeding"),
                       mutate(rsquared(sem_hab), sem = "habitat"))
model_fit <- imap_dfr(list(breeding = sem_bre, habitat = sem_hab), \(s, nm) {
  f <- fisherC(s)
  tibble(sem = nm, Fisher.C = f$Fisher.C, df = f$df, P.Value = f$P.Value)
})

# ---- diagnostics ----
lm_list <- list(food_breeding_PC1 = m_food_bre1, food_arrival_PC1 = m_food_hab1,
                food_arrival_PC2 = m_food_hab2)
lm_diagnostics <- imap_dfr(lm_list, \(m, nm) {
  sw <- shapiro.test(resid(m))
  bpt <- bptest(m)
  tibble(model = nm, n = nobs(m), shapiro_W = round(unname(sw$statistic), 4),
         shapiro_p = round(sw$p.value, 4), bp = round(unname(bpt$statistic), 3),
         bp_df = unname(bpt$parameter), bp_p = round(bpt$p.value, 4))
})

pearson <- sum(residuals(m_bre_resp, type = "pearson")^2)
poisson_diag <- tibble(model = "total_fledglings", converged = m_bre_resp$converged,
                       dispersion = round(pearson / df.residual(m_bre_resp), 3),
                       p_overdispersion = round(pchisq(pearson, df.residual(m_bre_resp),
                                                       lower.tail = FALSE), 4))

p_fit <- fitted(m_hab_resp)
n_pred <- length(coef(m_hab_resp)) - 1L
binomial_diag <- tibble(model = "breeding_bin", converged = m_hab_resp$converged,
                        fitted_min = round(min(p_fit), 4), fitted_max = round(max(p_fit), 4),
                        no_separation = min(p_fit) > 1e-6 & max(p_fit) < 1 - 1e-6,
                        events = sum(pred_hab$breeding_bin), predictors = n_pred,
                        events_per_predictor = round(sum(pred_hab$breeding_bin) / n_pred, 1))

sem_models <- list(
  bre_total_fledglings = list(m = m_bre_resp, ids = pred_bre$Site_ID),
  bre_food_breeding_PC1 = list(m = m_food_bre1, ids = pred_bre$Site_ID),
  hab_breeding_bin = list(m = m_hab_resp, ids = pred_hab$Site_ID),
  hab_food_arrival_PC1 = list(m = m_food_hab1, ids = pred_hab$Site_ID),
  hab_food_arrival_PC2 = list(m = m_food_hab2, ids = pred_hab$Site_ID)
)

influence_sem <- imap_dfr(sem_models, \(x, nm) {
  cd <- cooks.distance(x$m)
  rs <- rstudent(x$m)
  tibble(model = nm, n = length(cd),
         max_cooks_d = round(max(cd), 3), site_max_cooks_d = x$ids[which.max(cd)],
         max_abs_rstudent = round(max(abs(rs)), 2), site_max_rstudent = x$ids[which.max(abs(rs))],
         flagged = max(cd) > 1 | max(abs(rs)) > 3)
})

moran_sem <- imap_dfr(sem_models, \(x, nm) {
  r <- moran_resid(x$m, x$ids)
  tibble(model = nm, observed = round(r$observed, 3), expected = round(r$expected, 3),
         sd = round(r$sd, 3), p_value = round(r$p.value, 4))
})

# ---- bootstrap CI of the food availability lm ----
boot_coefs <- function(m, B) {
  dat <- model.frame(m)
  draws <- replicate(B, coef(lm(formula(m), data = dat[sample(nrow(dat), replace = TRUE), ])))
  tibble(term = rownames(draws), estimate = unname(coef(m)[rownames(draws)]),
         boot_lo = apply(draws, 1, quantile, probs = 0.025, na.rm = TRUE),
         boot_hi = apply(draws, 1, quantile, probs = 0.975, na.rm = TRUE))
}
set.seed(42)
boot_ci <- list(m_food_hab1 = m_food_hab1, m_food_hab2 = m_food_hab2) |>
  imap_dfr(\(m, nm) mutate(boot_coefs(m, n_boot), model = nm)) |>
  filter(term != "(Intercept)") |>
  mutate(ci_excludes_0 = boot_lo > 0 | boot_hi < 0,
         across(c(estimate, boot_lo, boot_hi), \(x) round(x, 3)))

# ---- Q-Q plot ----
qq_vars <- c("food_arrival_PC1", "food_arrival_PC2", "food_breeding_PC1")
qq_labels <- set_names(paste0(lbl(qq_vars), c(" (arrival)", " (arrival)", " (breeding)")), qq_vars)
qq_data <- imap_dfr(lm_list, \(m, nm) tibble(model = qq_labels[[nm]], resid = resid(m))) |>
  mutate(model = factor(model, levels = qq_labels))

p_qq <- ggplot(qq_data, aes(sample = resid)) +
  stat_qq(alpha = 0.7) +
  stat_qq_line(colour = col_food, linewidth = 0.8) +
  facet_wrap(~ model, ncol = 3) +
  xlab("theoretical quantiles") + ylab("sample quantiles") +
  theme_thesis(14)

ggsave(file.path(output_dir, "qqplots_sem_resid.png"), p_qq,
       width = 10, height = 4.5, device = agg_png, dpi = 300)

# ---- path diagrams ----
legend_txt <- paste0("Edge labels are standardised path coefficients. Solid lines, p &lt; ", alpha,
                     ". Dashed lines, p &#8805; ", alpha,
                     ". R&#178; of the GLMs is the Nagelkerke R&#178;.")

dot_nodes <- function(nodes) {
  r2 <- sprintf("R&#178; = %.2f", r_squared$R.squared[match(nodes$var, r_squared$Response)])
  style <- c(exo = ' fillcolor=white', med = paste0(' fillcolor="', fill_food, '"'),
             resp = paste0(' fillcolor="', col_resp, '" fontcolor=white'))
  sprintf(' %s [label="%s"%s]\n', nodes$id,
          ifelse(nodes$role == "exo", nodes$label, paste0(nodes$label, "\\n", r2)),
          style[nodes$role]) |>
    paste(collapse = "")
}

dot_edges <- function(cf, nodes, resp) {
  ids <- set_names(nodes$id, nodes$var)
  col <- ifelse(cf$Response == resp, col_resp, col_food)
  sig <- cf$P.Value < alpha
  sprintf(' %s -> %s [color="%s" fontcolor="%s" style=%s penwidth=%.1f label="%s"]\n',
          ids[cf$Predictor], ids[cf$Response], col, col,
          ifelse(sig, "solid", "dashed"), ifelse(sig, 2.4, 1.2),
          gsub("-", "&#8722;", sprintf("%.2f", cf$Std.Estimate))) |>
    paste(collapse = "")
}

save_path_diagram <- function(name, cf, nodes, resp, file, width) {
  blocks <- c("exo", "med", "resp") |>
    map_chr(\(r) paste0(" { rank=same\n", dot_nodes(filter(nodes, role == r)), " }\n"))
  dot <- paste0(
    "digraph ", name, " {\n",
    ' graph [rankdir=LR fontname="Calibri" bgcolor=white nodesep=0.45 ranksep=1.3\n',
    '  labelloc=b fontsize=11 label="', legend_txt, '"]\n',
    ' node [fontname="Calibri" fontsize=12 shape=box style=filled width=2.1 height=0.62]\n',
    ' edge [fontname="Calibri" fontsize=11]\n',
    paste(blocks, collapse = ""), dot_edges(cf, nodes, resp), "}\n")
  tmp <- tempfile(fileext = ".svg")
  grViz(dot) |> export_svg() |> writeLines(tmp)
  rsvg_png(tmp, file = file.path(output_dir, file), width = width)
}

add_labels <- function(nodes) {
  mutate(nodes, label = gsub("\u00b2", "&#178;", node_lab(var), fixed = TRUE), .before = role)
}

nodes_hab <- tribble(
  ~var, ~id, ~role,
  "site_structure_arrival_PC1", "ss", "exo",
  "land_use_arrival_PC1", "lu1", "exo",
  "land_use_arrival_PC2", "lu2", "exo",
  "veg_shannon", "vs", "exo",
  "n_perches", "np", "exo",
  "food_arrival_PC1", "fp1", "med",
  "food_arrival_PC2", "fp2", "med",
  "breeding_bin", "terr", "resp"
) |> add_labels()
nodes_bre <- tribble(
  ~var, ~id, ~role,
  "weather_PC1", "w1", "exo",
  "weather_PC2", "w2", "exo",
  "pesticides_PC1", "pe", "exo",
  "veg_shannon", "vs", "exo",
  "n_perches_25sqm", "np", "exo",
  "food_breeding_PC1", "fb", "med",
  "total_fledglings", "fl", "resp"
) |> add_labels()

save_path_diagram("hab_sem", coef_hab, nodes_hab, "breeding_bin", "sem_path_hab.png", 2200)
save_path_diagram("bre_sem", coef_bre, nodes_bre, "total_fledglings", "sem_path_bre.png", 1900)

# ---- export ----
write_xlsx(list(path_coeffs_breeding = coef_bre, path_coeffs_habitat = coef_hab,
                indirect_breeding = indirect_bre, indirect_habitat = indirect_hab,
                r_squared = as_tibble(r_squared), model_fit = model_fit,
                lm_diagnostics = lm_diagnostics, poisson_dispersion = poisson_diag,
                binomial_separation = binomial_diag, influence_sem = influence_sem,
                moran_sem = moran_sem, bootstrap_ci_hab = boot_ci),
           file.path(output_dir, "sem_results.xlsx"))
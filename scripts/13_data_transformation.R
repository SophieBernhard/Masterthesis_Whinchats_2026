# 13_data_transformation.R ----
# Check skewness (g1) of the raw predictors, compare candidate transformations, transform where required
# Flag rule: g1 > 2 OR (Shapiro-Wilk p < 0.05 AND g1 > 1), exceptions are listed in kept_raw
# All transformed columns carry the suffix _log (log, log1p or log(x - 1), see decision table)
# Requires: data_clean/master_habitat.rds and master_breeding.rds (script 11)
# Output: data_clean/master_habitat_tf.rds and master_breeding_tf.rds (read by script 14)

# ---- packages ----
library("tidyverse")
library("writexl")
library("here")

# ---- paths ----
data_dir <- here("data_clean")
output_dir <- here("output")

# ---- settings ----
skew_hi <- 2 # g1 above this value flags a variable on skewness alone
skew_vis <- 1 # g1 above this value flags a variable together with a significant Shapiro-Wilk test
alpha <- 0.05 # significance level of the Shapiro-Wilk test
cand_thr <- 0.5 # |g1| above this value enters the comparison of candidate transformations
kept_raw <- c("pct_arable", "pct_intensive", "n_perches") # flagged but kept raw (zero-inflated land use, single predictor)
ds <- c(hab = "n = 40", bre = "n = 19") # dataset label used in all tables

# ---- load ----
hab <- readRDS(file.path(data_dir, "master_habitat.rds"))
bre <- readRDS(file.path(data_dir, "master_breeding.rds"))
dat <- list(hab = hab, bre = bre)

# ---- input checks ----
stopifnot(nrow(hab) == 40, ncol(hab) == 25, !anyDuplicated(hab$Site_ID))
stopifnot(nrow(bre) == 19, ncol(bre) == 36, !anyDuplicated(bre$Site_ID))

# ---- skewness check (raw) ----
skew_g1 <- function(x) {
  x <- x[!is.na(x)]
  z <- x - mean(x)
  mean(z^3) / mean(z^2)^1.5
}
skew_table <- function(d, vars, nm) {
  d |>
    select(all_of(vars)) |>
    pivot_longer(everything()) |>
    summarise(n = sum(!is.na(value)),
              n_zero = sum(value == 0, na.rm = TRUE),
              skew = round(skew_g1(value), 2),
              sw_p = if (n_distinct(value, na.rm = TRUE) >= 3) round(shapiro.test(value)$p.value, 4) else NA_real_,
              .by = name) |>
    mutate(dataset = nm)
}
not_tested <- c("Site_ID", "breeding_bin", "total_fledglings", "n_pairs_total", "fledglings_na")
vars_hab <- hab |> select(where(is.numeric), -any_of(not_tested)) |> names()
vars_bre <- bre |> select(where(is.numeric), -any_of(not_tested)) |> names()
skew_raw <- bind_rows(skew_table(hab, vars_hab, ds[["hab"]]), skew_table(bre, vars_bre, ds[["bre"]]))
print(skew_raw, n = Inf)

# ---- comparison of transformations against raw ----
fns <- list(raw = \(x) x,
            log = \(x) if (all(x > 0, na.rm = TRUE)) log(x) else x * NA,
            log1p = \(x) if (all(x >= 0, na.rm = TRUE)) log1p(x) else x * NA,
            sqrt = \(x) if (all(x >= 0, na.rm = TRUE)) sqrt(x) else x * NA,
            logm1 = \(x) if (all(x > 1, na.rm = TRUE)) log(x - 1) else x * NA)
cmp_tf <- function(d, nm) {
  v <- skew_raw |> filter(dataset == nm, abs(skew) > cand_thr, !str_detect(name, "nmds")) |> pull(name)
  map_dfr(names(fns), \(f) d |> select(all_of(v)) |> mutate(across(everything(), fns[[f]])) |>
            skew_table(v, nm) |> mutate(fun = f))
}
skew_cand <- bind_rows(cmp_tf(hab, ds[["hab"]]), cmp_tf(bre, ds[["bre"]])) |>
  filter(fun != "logm1" | name == "rumple") |>
  mutate(gain = round(abs(skew[fun == "raw"]) - abs(skew), 2), .by = c(dataset, name)) |>
  select(dataset, name, fun, skew, sw_p, gain) |>
  arrange(dataset, name, factor(fun, names(fns)))
print(skew_cand, n = Inf)

# ---- transformation decision ----
tf <- tribble(
  ~dataset, ~variable, ~fun, ~reason,
  "hab", "mean_h", "log", "g1 > 2",
  "hab", "terrain_sd", "log", "g1 > 2",
  "hab", "pct_settlement", "log1p", "g1 > 2",
  "hab", "pct_seminat", "log1p", "Shapiro-Wilk significant, right-skewed",
  "hab", "rumple", "logm1", "Shapiro-Wilk significant, right-skewed",
  "bre", "mean_h", "log", "Shapiro-Wilk significant, right-skewed",
  "bre", "pct_settlement", "log1p", "g1 > 2",
  "bre", "pct_seminat", "log1p", "Shapiro-Wilk significant, right-skewed",
  "bre", "biomass_mean_breeding", "log", "Shapiro-Wilk significant, right-skewed"
) |>
  mutate(new = paste0(variable, "_log"))
tf_fun <- c(log = log, log1p = log1p, logm1 = \(x) log(x - 1))
print(tf, n = Inf)

# ---- decision checks ----
# domain of each chosen function
dom_ok <- c(log = \(x) all(x > 0), log1p = \(x) all(x >= 0), logm1 = \(x) all(x > 1))
for (i in seq_len(nrow(tf))) {
  x <- na.omit(dat[[tf$dataset[i]]][[tf$variable[i]]])
  stopifnot(is.numeric(x), dom_ok[[tf$fun[i]]](x))
}
# every flagged variable is either transformed or listed in kept_raw
flagged <- skew_raw |>
  filter(skew > skew_hi | (coalesce(sw_p, 1) < alpha & skew > skew_vis)) |>
  mutate(key = names(ds)[match(dataset, ds)],
         transformed = paste(key, name) %in% paste(tf$dataset, tf$variable))
print(flagged, n = Inf)
stopifnot(all(paste(tf$dataset, tf$variable) %in% paste(flagged$key, flagged$name)),
          setequal(flagged$name[!flagged$transformed], kept_raw))
# every chosen function has the smallest |g1| among the valid candidates
best <- skew_cand |>
  filter(!is.nan(skew)) |>
  summarise(best = min(abs(skew)), .by = c(dataset, name))
chosen <- tf |>
  mutate(dataset = unname(ds[dataset])) |>
  left_join(skew_cand, by = c("dataset", "variable" = "name", "fun")) |>
  left_join(best, by = c("dataset", "variable" = "name"))
stopifnot(!anyNA(chosen$skew), all(abs(chosen$skew) <= chosen$best))

# ---- transformation ----
transform_data <- function(d, nm) {
  for (i in which(tf$dataset == nm)) {
    res <- tf_fun[[tf$fun[i]]](d[[tf$variable[i]]])
    stopifnot(all(is.finite(res[!is.na(d[[tf$variable[i]]])])))
    d[[tf$new[i]]] <- res
    d[[tf$variable[i]]] <- NULL
  }
  d
}
hab_tf <- transform_data(hab, "hab")
bre_tf <- transform_data(bre, "bre")

# ---- result checks ----
check_tf <- function(d, d_tf, nm) {
  v <- tf$variable[tf$dataset == nm]
  new <- tf$new[tf$dataset == nm]
  keep <- setdiff(names(d), v)
  stopifnot(ncol(d_tf) == ncol(d), identical(d_tf$Site_ID, d$Site_ID),
            !any(v %in% names(d_tf)), all(new %in% names(d_tf)),
            all(unname(colSums(is.na(d_tf[new]))) == unname(colSums(is.na(d[v])))),
            isTRUE(all.equal(d_tf[keep], d[keep])))
}
check_tf(hab, hab_tf, "hab")
check_tf(bre, bre_tf, "bre")

# ---- skewness check (after transformation, next to raw) ----
post <- \(d, key) {
  skew_table(d, tf$new[tf$dataset == key], ds[[key]]) |>
    mutate(variable = tf$variable[tf$dataset == key][match(name, tf$new[tf$dataset == key])]) |>
    left_join(skew_raw |> filter(dataset == ds[[key]]) |> select(variable = name, skew_raw = skew, sw_raw = sw_p),
              by = "variable")
}
skew_post <- bind_rows(post(hab_tf, "hab"), post(bre_tf, "bre"))
print(skew_post, n = Inf)
stopifnot(nrow(skew_post) == nrow(tf), all(abs(skew_post$skew) < abs(skew_post$skew_raw)))

# ---- summary ----
cat("Variables tested: ", nrow(skew_raw), " (", sum(skew_raw$dataset == ds[["hab"]]), " in n = 40, ",
    sum(skew_raw$dataset == ds[["bre"]]), " in n = 19)\n", sep = "")
cat("Flagged: ", nrow(flagged), ", transformed: ", nrow(tf), ", kept raw: ", length(kept_raw), " (",
    paste(kept_raw, collapse = ", "), ")\n", sep = "")
cat("master_habitat_tf: ", nrow(hab_tf), " x ", ncol(hab_tf), ", master_breeding_tf: ",
    nrow(bre_tf), " x ", ncol(bre_tf), "\n", sep = "")

# ---- save ----
saveRDS(hab_tf, file.path(data_dir, "master_habitat_tf.rds"))
saveRDS(bre_tf, file.path(data_dir, "master_breeding_tf.rds"))
 write_xlsx(hab_tf, file.path(data_dir, "master_habitat_tf.xlsx"))
 write_xlsx(bre_tf, file.path(data_dir, "master_breeding_tf.xlsx"))
 write_xlsx(list(raw = skew_raw, candidates = skew_cand, transformed = skew_post, decision = tf),
            file.path(output_dir, "13_skewness_transformation.xlsx"))
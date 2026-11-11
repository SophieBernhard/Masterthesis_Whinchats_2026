# 01_food_data_cleaning.R ----
# Load and clean COI metabarcoding data, then compute food diversity variables per site and phase

# ---- packages ----
library("readxl")
library("tidyverse")
library("writexl")
library("here")
library("vegan")

# ---- paths ----
food_data_raw <- here("data_raw/Mueller_EC24-314_COI_results.xlsx")
output_dir <- here("data_clean")

# ---- settings ----
min_identity <- 0.98 # BOLD hit and grade identity (proportion)
min_reads <- 1       # reads per site and OTU needed to count as present
k_chosen <- 4        # NMDS dimensions, see stress_by_k in the diagnostics output
min_sites <- 10      # families tested against NMDS axes need at least this many sites

# ---- metadata (columns A:Z) ----
food_meta_raw <- read_excel(
  food_data_raw, skip = 3, guess_max = 6254, .name_repair = "minimal"
) |>
  select(
    OTU_ID_raw = "OTU_ID;cluster_size",
    BOLD_HIT_ID = "BOLD_HIT%ID",
    BOLD_Grade_ID = "BOLD_Grade%ID",
    BOLD_Phylum, BOLD_Class, BOLD_Order, BOLD_Family, BOLD_Genus, BOLD_Species
  ) |>
  mutate(OTU_ID = str_extract(OTU_ID_raw, "^[^;]+")) |>
  select(-OTU_ID_raw)

# ---- sample names (row 3, columns AA:FZ) ----
food_sample_names <- read_excel(
  food_data_raw, range = "AA3:FZ3", col_names = FALSE, .name_repair = "minimal"
) |>
  unlist() |>
  as.character()

# ---- read count matrix (row 5 onwards, columns AA:FZ) ----
food_reads_raw <- read_excel(
  food_data_raw,
  range = paste0("AA5:FZ", 4 + nrow(food_meta_raw)),
  col_names = FALSE, .name_repair = "minimal"
) |>
  as.matrix()
colnames(food_reads_raw) <- food_sample_names

# ---- filter: high identity BOLD hits, Insecta and Arachnida only ----
keep <- which(
  food_meta_raw$BOLD_HIT_ID >= min_identity &
    food_meta_raw$BOLD_Grade_ID >= min_identity &
    food_meta_raw$BOLD_Class %in% c("Insecta", "Arachnida")
)
food_filtered <- food_meta_raw[keep, ]
food_reads_filtered <- food_reads_raw[keep, ]
stopifnot(!anyNA(food_filtered$BOLD_Family), !anyNA(food_reads_filtered))

# ---- sample metadata (site ID and month) ----
site_info <- tibble(sample = food_sample_names) |>
  mutate(
    Site_ID = str_extract(sample, "^[^_]+"),
    month = month(as.Date(str_extract(sample, "\\d{2}_\\d{2}_\\d{2}$"),
                          format = "%d_%m_%y"))
  )

# ---- presence matrices per phase ----
# Sites in rows, OTUs or families in columns
make_pa <- function(phase_months) {
  smp <- site_info |> filter(month %in% phase_months)
  site_reads <- rowsum(t(food_reads_filtered[, smp$sample, drop = FALSE]), smp$Site_ID)
  otu_pa <- (site_reads >= min_reads) * 1L
  fam_pa <- (t(rowsum(t(otu_pa), food_filtered$BOLD_Family)) > 0) * 1L
  list(otu = otu_pa, fam = fam_pa)
}

# ---- NMDS diagnostics ----
nmds_diagnostics <- function(phase_label, dist_fam, nmds_fit) {
  region <- str_extract(attr(dist_fam, "Labels"), "^[A-Z]+")
  set.seed(42)
  bd_p <- permutest(betadisper(dist_fam, region), permutations = 999)$tab[1, "Pr(>F)"]
  stress_k <- map_dbl(1:6, \(k) {
    set.seed(42)
    metaMDS(dist_fam, k = k, trymax = 100, halfchange = TRUE, trace = FALSE)$stress
  })
  list(
    stress_by_k = tibble(phase = phase_label, k = 1:6, stress = stress_k),
    fit = tibble(
      phase = phase_label, k = nmds_fit$ndim, stress = nmds_fit$stress,
      converged_runs = as.integer(nmds_fit$converged),
      nonmetric_r2 = 1 - nmds_fit$stress^2,
      n_sites = attr(dist_fam, "Size"), n_dissimilarities = length(dist_fam),
      n_parameters = nmds_fit$ndim * attr(dist_fam, "Size"),
      betadisper_p_region = bd_p
    )
  )
}

# ---- food diversity function ----
compute_food_diversity <- function(phase_label, phase_months) {
  pa <- make_pa(phase_months)
  stopifnot(all(rowSums(pa$fam) > 0))
  sites <- rownames(pa$fam)
  
  dist_fam <- vegdist(pa$fam, method = "jaccard", binary = TRUE)
  # halfchange = TRUE fixes the axis scaling across vegan versions
  set.seed(42)
  nmds_fit <- metaMDS(dist_fam, k = k_chosen, trymax = 100, halfchange = TRUE, trace = FALSE)
  stopifnot(as.integer(nmds_fit$converged) > 0)
  message(phase_label, " NMDS stress (k=", k_chosen, "): ", round(nmds_fit$stress, 4))
  nmds_scores <- scores(nmds_fit, display = "sites")
  
  # Spearman correlation of each family with each NMDS axis, FDR within phase and axis
  fam_corr <- map_dfr(colnames(pa$fam), \(fam) {
    n_present <- sum(pa$fam[, fam])
    if (n_present < min_sites || n_present == length(sites)) return(NULL)
    map_dfr(seq_len(k_chosen), \(ax) {
      tst <- cor.test(pa$fam[, fam], nmds_scores[sites, ax], method = "spearman", exact = FALSE)
      tibble(phase = phase_label, family = fam, axis = paste0("NMDS", ax),
             n_sites = length(sites), n_present = n_present,
             rho = unname(tst$estimate), p_raw = tst$p.value)
    })
  }) |>
    group_by(axis) |>
    mutate(p_fdr = p.adjust(p_raw, method = "fdr")) |>
    ungroup()
  
  # Mean Jaccard distance of each site to all other sites
  dm <- as.matrix(dist_fam)
  diag(dm) <- NA
  
  list(
    diversity = tibble(
      Site_ID = sites,
      !!paste0("otu_richness_", phase_label) := rowSums(pa$otu)[sites],
      !!paste0("food_nmds1_", phase_label) := nmds_scores[sites, 1],
      !!paste0("food_nmds2_", phase_label) := nmds_scores[sites, 2],
      !!paste0("food_nmds3_", phase_label) := nmds_scores[sites, 3],
      !!paste0("food_nmds4_", phase_label) := nmds_scores[sites, 4],
      !!paste0("food_beta_mean_", phase_label) := rowMeans(dm, na.rm = TRUE)[sites]
    ),
    fam_corr = fam_corr,
    diagnostics = nmds_diagnostics(phase_label, dist_fam, nmds_fit)
  )
}

# ---- run for arrival (May) and breeding (June and July) ----
div_arrival <- compute_food_diversity("arrival", phase_months = 5)
div_breeding <- compute_food_diversity("breeding", phase_months = c(6, 7))

# ---- save ----
food_diversity_split <- left_join(div_arrival$diversity, div_breeding$diversity, by = "Site_ID")
saveRDS(food_diversity_split, file.path(output_dir, "food_diversity_split.rds"))

fam_corr_all <- bind_rows(div_arrival$fam_corr, div_breeding$fam_corr)
 write_xlsx(
  list(all_tests = fam_corr_all, significant = filter(fam_corr_all, p_fdr < 0.05)),
  file.path(output_dir, "food_nmds_family_correlations.xlsx")
 )

 write_xlsx(
  list(
    stress_by_k = bind_rows(div_arrival$diagnostics$stress_by_k, div_breeding$diagnostics$stress_by_k),
    final_fit = bind_rows(div_arrival$diagnostics$fit, div_breeding$diagnostics$fit)
  ),
  file.path(output_dir, "food_nmds_diagnostics.xlsx")
 )
 write_xlsx(food_diversity_split, file.path(output_dir, "food_diversity_split.xlsx"))
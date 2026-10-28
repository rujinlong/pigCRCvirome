# Compare regenerated results with the source-data workbook submitted with the Letter.
here::i_am("scripts/verify.R")
suppressPackageStartupMessages(library(dplyr))
input <- function(file) here::here("analyses", "data", "00-raw", file)
read_input <- function(name) readr::read_tsv(input(paste0(name, ".tsv.gz")), show_col_types = FALSE)
community <- readRDS(here::here("analyses", "data", "01-community", "community_results.rds"))
life <- readRDS(here::here("analyses", "data", "02-lifestyle", "lifestyle_results.rds"))
host <- readRDS(here::here("analyses", "data", "03-host-associations", "host_results.rds"))
month <- c(Aug = 3, Sep = 4, Dec = 7)
site <- c(fecal = "feces", cecum = "cecal contents", colon = "colonic contents")
sample_type <- c(fecal = "feces", digesta = "gut contents")
model <- c(standard = "primary", depth_adjusted = "mapped-depth adjusted")
test <- c(fisher = "Fisher's exact test", quality_VC = "logistic model adjusted for length and completeness")
kw_columns <- function(df) transmute(df, n_pigs = n, kruskal_wallis_h = statistic,
  p_chi_square = p_chisq, p_permutation = p_perm, n_permutations = nperm)

tables <- list(
  D1_libraries = read_input("samples") %>% transmute(library_id = id_virome, pig_id = host_id,
    diet = Diet, study_month = unname(month[Month_raw]), site = unname(site[Location])),
  D2_ordination = community$ordination %>% transmute(pig_id = host_id, diet = Diet,
    pcoa1 = PCoA1, pcoa2 = PCoA2, pcoa1_variance_percent = axis1_percent,
    pcoa2_variance_percent = axis2_percent),
  D3_community_tests = community$community %>% transmute(study_month = unname(month[Month_raw]),
    site = unname(site[Location]), n_pigs = n, n_votus = n_votu, r_squared, pseudo_f,
    p_permutation = p_perm, permdisp_p = dispersion_p, n_permutations = nperm),
  D4_paired_change = community$paired_change %>% transmute(n_pigs = n, r_squared, pseudo_f, p_permutation = p_perm),
  D5_litter_tests = community$litter %>% transmute(
    comparison = recode(analysis, endpoint_fecal = "endpoint feces", paired_change = "month-7 minus month-4 change"),
    litter_blocks = recode(scenario, known_litter_only = "pigs with known litter",
      unknown_singletons = "unknown-litter pigs as single-pig blocks"),
    n_pigs = n, n_blocks, n_label_assignments = n_assignments,
    n_assignments_f_at_least_observed = tail_count, p_minimum = p_min, r_squared, pseudo_f, p_exact),
  D6_lifestyle_enrichment = life$enrichment %>% transmute(sample_type = unname(sample_type[compartment]),
    model = unname(.env$model[model]), test = unname(test[inference]), n_tested_votus = n_tested,
    n_diet_associated = n_da, diet_associated_temperate = da_temperate, diet_associated_other = da_other,
    other_tested_temperate = non_da_temperate, other_tested_other = non_da_other,
    diet_associated_temperate_percent = da_percent, other_tested_temperate_percent = non_da_percent,
    odds_ratio, ci95_low = ci_low, ci95_high = ci_high, p, n_vcontact3_clusters = n_clusters, q_bh = q_compartments),
  D7_temperate_share_test = kw_columns(life$share_test),
  D8_bai_host_counts = host$census %>% transmute(host_family = top_genome_family,
    host_genus = top_genome_genus, predicted_lifestyle = lifestyle, n_votus = n_votu),
  D9_candidate_detection = host$detection %>% transmute(pig_id = host_id, diet = Diet,
    study_month = unname(month[Month_raw]), vOTU_id, tpm, covered_fraction, detected),
  D10_candidate_genes = host$genes %>% transmute(vOTU_id, gene_id = protein_id,
    start = start_position, end = end_position, strand = ifelse(strandedness == 1, "forward", "reverse"),
    functional_class, highlighted_annotation = label),
  D11_endpoint_pigs = host$endpoint %>% transmute(pig_id = host_id, diet = Diet,
    baseline_distal_colon_polyps = pre_n_40cm, progressing_polyps = Growing_Polyps,
    desulfovibrionaceae_phage_percent = 100 * desulfo_share),
  D12_diet_medians = host$medians %>% transmute(diet = Diet, n_pigs = n,
    median_desulfovibrionaceae_phage_percent = desulfo_median_percent),
  D13_diet_test = kw_columns(host$diet_test),
  D14_pairwise_tests = bind_cols(select(host$pairwise, contrast), kw_columns(host$pairwise), select(host$pairwise, q_bh)),
  D15_conditional_tests = host$conditional %>% transmute(covariate, n_pigs = nrow(host$endpoint),
    statistic_name = "F", statistic = f, df_numerator = df_diet, residual_df = df_residual,
    p_permutation = p_perm, n_permutations = nperm,
    permutation = "Unrestricted Freedman-Lane residual permutation"),
  D16_desulfo_genera = host$genera %>% transmute(host_genus, n_votus = n_votu),
  D17_family_screen = host$family_screen %>% filter(host_family == "Desulfovibrionaceae") %>%
    transmute(host_family, n_pigs = n, n_covariates, residual_df, partial_rho = rho,
      p_permutation = p_perm, ci95_low = ci_low, ci95_high = ci_high,
      n_permutations = nperm, q_bh = q_bh_perm),
  D18_origin_test = kw_columns(host$origin_test),
  D19_origin_medians = host$origin_medians %>% transmute(diet = Diet, n_pigs = n,
    median_desulfovibrionaceae_phage_percent = median_percent,
    median_detected_votus = median_detected_cf20, min_detected_votus = min_detected_cf20,
    max_detected_votus = max_detected_cf20)
)
keys <- list("library_id", "pig_id", c("study_month", "site"), "n_pigs", "comparison",
  c("sample_type", "model", "test"), "n_pigs", c("host_family", "host_genus", "predicted_lifestyle"),
  c("vOTU_id", "pig_id", "study_month"), "gene_id", "pig_id", "diet", "n_pigs", "contrast",
  "covariate", "host_genus", "host_family", "n_pigs", "diet")
sort_table <- function(df, key) {
  df <- as.data.frame(lapply(df, function(x) if (is.factor(x)) as.character(x) else x))
  df <- df[do.call(order, df[key]), , drop = FALSE]
  rownames(df) <- NULL
  df
}
failures <- character()
for (i in seq_along(tables)) {
  name <- names(tables)[i]
  actual <- sort_table(tables[[i]], keys[[i]])
  expected <- sort_table(readxl::read_xlsx(input("submitted_source_data.xlsx"), sheet = name), keys[[i]])
  if (!identical(names(actual), names(expected)) || nrow(actual) != nrow(expected)) {
    failures <- c(failures, paste(name, "table dimensions or column names differ"))
    next
  }
  # Eigenvector signs are arbitrary; compare ordination axes up to reflection.
  if (name == "D2_ordination") for (axis in c("pcoa1", "pcoa2")) {
    if (sum(actual[[axis]] * expected[[axis]]) < 0) actual[[axis]] <- -actual[[axis]]
  }
  differences <- vapply(names(expected), function(column) {
    x <- actual[[column]]; y <- expected[[column]]
    if (!identical(is.na(x), is.na(y))) return(TRUE)
    keep <- !is.na(y)
    if (is.numeric(y)) any(abs(x[keep] - y[keep]) > 1e-9 * pmax(1, abs(y[keep]))) else
      !identical(as.character(x[keep]), as.character(y[keep]))
  }, logical(1))
  if (any(differences)) failures <- c(failures, paste(name, paste(names(differences)[differences], collapse = ", ")))
  else message("PASS ", name, " (", nrow(actual), " rows)")
}

# These rounded values are reported in the submitted supplementary methods.
pairwise <- community$pairwise
stopifnot(identical(pairwise$contrast, c("CTRL vs RL", "CTRL vs COL", "RL vs COL")),
          identical(round(pairwise$r_squared, 3), c(.214, .166, .131)),
          identical(round(pairwise$p_perm, 3), c(.009, .003, .057)),
          identical(round(pairwise$q_bh, 3), c(.013, .008, .057)))
message("PASS supplementary community pairwise contrasts (3 rows)")

manifest <- readr::read_tsv(input("checksums.tsv"), show_col_types = FALSE)
hashes <- vapply(manifest$file, function(file) digest::digest(input(file), algo = "sha256", file = TRUE), character(1))
stopifnot(identical(unname(hashes), manifest$sha256))
for (figure in c("Figure_1", "Figure_2")) for (format in c("png", "svg", "pdf", "tiff")) {
  path <- here::here("analyses", "data", "04-figures", paste0(figure, ".", format))
  stopifnot(file.exists(path), file.info(path)$size > 0)
}
if (length(failures)) stop("Submitted-data comparison failed:\n", paste(failures, collapse = "\n"))
message("Verified all 19 submitted source-data tables, input checksums, and both figure exports.")

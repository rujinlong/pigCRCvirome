# Statistical tests used in the Gut Letter.

make_relative_abundance <- function(mat_abundance) {
  stopifnot(all(is.finite(mat_abundance)), all(mat_abundance >= 0),
            all(colSums(mat_abundance) > 0))
  sweep(mat_abundance, 2, colSums(mat_abundance), "/")
}

make_permutation_indices <- function(n, nperm, strata = NULL) {
  if (is.null(strata)) return(replicate(nperm, sample.int(n)))
  stopifnot(length(strata) == n, !anyNA(strata))
  groups <- split(seq_len(n), strata)
  replicate(nperm, {
    idx <- seq_len(n)
    for (group in groups) idx[group] <- group[sample.int(length(group))]
    idx
  })
}

test_rank_association <- function(x, y, df_covariates = NULL,
                                  strata = NULL, nperm = 9999L) {
  keep <- is.finite(x) & is.finite(y)
  if (!is.null(df_covariates)) keep <- keep & stats::complete.cases(df_covariates)
  if (!is.null(strata)) keep <- keep & !is.na(strata)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  if (is.null(df_covariates)) {
    mat_covariates <- matrix(1, nrow = n, ncol = 1)
  } else {
    df_covariates <- droplevels(as.data.frame(df_covariates[keep, , drop = FALSE]))
    df_covariates[] <- lapply(df_covariates, function(v)
      if (is.numeric(v)) rank(v, ties.method = "average") else v)
    mat_covariates <- stats::model.matrix(~ ., df_covariates)
  }
  qr_covariates <- qr(mat_covariates)
  n_covariates <- qr_covariates$rank - 1L
  residual_df <- n - n_covariates - 2L
  if (n < 4L || residual_df <= 0L || length(unique(x)) < 2L ||
      length(unique(y)) < 2L) {
    return(data.frame(n, n_covariates, residual_df, rho = NA_real_,
                      p_t = NA_real_, p_perm = NA_real_,
                      ci_low = NA_real_, ci_high = NA_real_, nperm = 0L))
  }
  rx <- qr.resid(qr_covariates, rank(x, ties.method = "average"))
  ry <- qr.resid(qr_covariates, rank(y, ties.method = "average"))
  rho <- sum(rx * ry) / sqrt(sum(rx^2) * sum(ry^2))
  rho <- min(max(rho, -1), 1)
  statistic <- rho * sqrt(residual_df / (1 - rho^2))
  p_t <- 2 * stats::pt(-abs(statistic), df = residual_df)
  p_perm <- NA_real_
  if (nperm > 0L) {
    mat_indices <- make_permutation_indices(n, nperm,
                                            if (is.null(strata)) NULL else strata[keep])
    # Freedman-Lane refits nuisance terms after each residual permutation.
    mat_permuted <- matrix(ry[mat_indices], nrow = n)
    mat_permuted <- qr.resid(qr_covariates, mat_permuted)
    rho_permuted <- as.vector(crossprod(rx, mat_permuted)) /
      sqrt(sum(rx^2) * colSums(mat_permuted^2))
    p_perm <- (1 + sum(abs(rho_permuted) >= abs(rho) - 1e-12)) / (nperm + 1)
  }
  z <- atanh(min(max(rho, -1 + 1e-12), 1 - 1e-12))
  se <- if (n - n_covariates - 3L > 0L) 1 / sqrt(n - n_covariates - 3L) else NA_real_
  data.frame(n, n_covariates, residual_df, rho, p_t, p_perm,
             ci_low = tanh(z - 1.96 * se), ci_high = tanh(z + 1.96 * se), nperm)
}

test_rank_groups <- function(x, groups, nperm = 9999L) {
  keep <- is.finite(x) & !is.na(groups)
  x <- x[keep]
  groups <- droplevels(factor(groups[keep]))
  ranks <- rank(x)
  n <- length(x)
  observed <- stats::kruskal.test(x, groups)
  mat_group <- stats::model.matrix(~ groups - 1)
  n_group <- colSums(mat_group)
  tie_sizes <- as.numeric(table(x))
  tie_correction <- 1 - sum(tie_sizes^3 - tie_sizes) / (n^3 - n)
  mat_indices <- make_permutation_indices(n, nperm)
  mat_sums <- crossprod(mat_group, matrix(ranks[mat_indices], nrow = n))
  permuted <- (12 / (n * (n + 1)) *
    colSums(sweep(mat_sums^2, 1, n_group, "/")) - 3 * (n + 1)) / tie_correction
  data.frame(n, statistic = unname(observed$statistic), p_chisq = observed$p.value,
             p_perm = (1 + sum(permuted >= observed$statistic - 1e-12)) / (nperm + 1),
             nperm)
}

test_community_groups <- function(mat_abundance, df_metadata,
                                   distance = "bray", nperm = 9999L) {
  stopifnot(ncol(mat_abundance) == nrow(df_metadata),
            !anyDuplicated(df_metadata$host_id))
  mat_relative <- t(make_relative_abundance(mat_abundance))
  dissimilarity <- vegan::vegdist(mat_relative, method = distance)
  fit <- vegan::adonis2(dissimilarity ~ Diet, data = df_metadata,
                        permutations = nperm)
  dispersion <- vegan::betadisper(dissimilarity, factor(df_metadata$Diet),
                                  type = "median", bias.adjust = TRUE)
  dispersion_test <- vegan::permutest(dispersion, permutations = nperm)
  data.frame(n = nrow(df_metadata), n_votu = nrow(mat_abundance),
             r_squared = fit$R2[1], pseudo_f = fit$F[1],
             p_perm = fit$`Pr(>F)`[1],
             dispersion_p = dispersion_test$tab$`Pr(>F)`[1], nperm)
}

enumerate_litter_labels <- function(diet, litter) {
  groups <- split(seq_along(diet), litter)
  choices <- lapply(groups, function(group) {
    permutations <- if (length(group) > 1L) {
      rbind(seq_along(group), as.matrix(permute::allPerms(length(group))))
    } else matrix(1L, 1, 1)
    unique(matrix(as.character(diet[group])[permutations],
                  nrow = nrow(permutations)))
  })
  df_combinations <- expand.grid(lapply(choices, function(x) seq_len(nrow(x))))
  mat_labels <- matrix(as.character(diet), nrow = nrow(df_combinations),
                        ncol = length(diet), byrow = TRUE)
  for (j in seq_along(groups)) {
    mat_labels[, groups[[j]]] <- choices[[j]][df_combinations[[j]], , drop = FALSE]
  }
  mat_labels
}

test_diet_adjusted <- function(y, groups, covariate, nperm = 19999L, seed = 138L) {
  y <- rank(y)
  mat_z <- cbind(1, rank(covariate))
  mat_x <- cbind(mat_z, stats::model.matrix(~ groups)[, -1, drop = FALSE])
  qr_z <- qr(mat_z); qr_x <- qr(mat_x)
  df_full <- length(y) - qr_x$rank
  df_diet <- qr_x$rank - qr_z$rank
  f_stat <- function(v) {
    rss_reduced <- colSums(qr.resid(qr_z, v)^2)
    rss_full <- colSums(qr.resid(qr_x, v)^2)
    ((rss_reduced - rss_full) / df_diet) / (rss_full / df_full)
  }
  residual <- qr.resid(qr_z, y)
  f_obs <- f_stat(matrix(y))
  set.seed(seed)
  idx <- replicate(nperm, sample.int(length(y)))
  f_perm <- f_stat((y - residual) + matrix(residual[idx], nrow = length(y)))
  tibble::tibble(f = f_obs, df_diet = df_diet, df_residual = df_full,
                 p_perm = (1 + sum(f_perm >= f_obs * (1 - 1e-12))) / (nperm + 1), nperm = nperm)
}

use_permutation_state <- function(states, name) {
  assign(".Random.seed", states[[name]], envir = .GlobalEnv)
  invisible(NULL)
}

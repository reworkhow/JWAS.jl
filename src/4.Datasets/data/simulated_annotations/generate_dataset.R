suppressPackageStartupMessages(library(utils))

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath("src/4.Datasets/data/simulated_annotations/generate_dataset.R", mustWork = TRUE)
}

dataset_dir <- dirname(script_path)

simulation_seed <- 123L
h2 <- 0.5
ncv <- 10L
S <- -1
n_gwas <- 400L
mt_n_shared <- 8L
mt_n_trait1_only <- 6L
mt_n_trait2_only <- 6L
mt_shared_rho <- 0.6
mt_residual_rho <- 0.3

# The packaged fixture is self-contained: regeneration reads the raw genotype
# source stored alongside this script rather than reaching back into temp/.
geno_path <- file.path(dataset_dir, "raw_genotypes.txt")
geno <- as.matrix(read.table(geno_path))
freq <- colMeans(geno) / 2
X <- geno[, freq > 0.01 & freq < 0.99, drop = FALSE]
nsnp <- ncol(X)

set.seed(simulation_seed)
cv <- sample(seq_len(nsnp), ncv)
freq_cv <- colMeans(X[, cv, drop = FALSE]) / 2
beta_cv <- rnorm(ncv, mean = 0, sd = sqrt((2 * freq_cv * (1 - freq_cv))^S))
g <- as.vector(X[, cv, drop = FALSE] %*% beta_cv)
e <- rnorm(nrow(X), 0, sqrt(var(g) * (1 - h2) / h2))
y <- g + e

X_t <- X[1:n_gwas, , drop = FALSE]
y_t <- y[1:n_gwas]

functional <- integer(nsnp)
random_anno <- integer(nsnp)
functional[cv] <- 1
noncausal <- setdiff(seq_len(nsnp), cv)
random_anno[sample(noncausal, 10)] <- 1

true_beta <- numeric(nsnp)
true_beta[cv] <- beta_cv

marker_ids <- paste0("m", seq_len(nsnp))
ids <- paste0("id_", seq_len(n_gwas))

geno_df <- data.frame(ID = ids, X_t, check.names = FALSE)
names(geno_df)[-1] <- marker_ids
write.csv(geno_df, file.path(dataset_dir, "genotypes.csv"), row.names = FALSE)
write.csv(data.frame(ID = ids, y1 = y_t), file.path(dataset_dir, "phenotypes.csv"), row.names = FALSE)
write.csv(
  data.frame(marker_id = marker_ids, functional = functional, random_anno = random_anno),
  file.path(dataset_dir, "annotations.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(marker_id = marker_ids, is_causal = seq_len(nsnp) %in% cv, true_effect = true_beta),
  file.path(dataset_dir, "truth.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# Native 2-trait extension for the packaged fixture.
#
# The 2-trait files are saved alongside the original single-trait files so the
# old fixture remains stable while multi-trait benchmarks can resolve a native
# two-trait dataset from the same packaged directory.
# ---------------------------------------------------------------------------
set.seed(simulation_seed + 1L)
mt_total_active <- mt_n_shared + mt_n_trait1_only + mt_n_trait2_only
if (mt_total_active >= nsnp) {
  stop("The requested number of 2-trait active markers must be smaller than the filtered marker count.")
}

mt_active <- sample(seq_len(nsnp), mt_total_active)
mt_shared <- mt_active[seq_len(mt_n_shared)]
mt_trait1_only <- mt_active[mt_n_shared + seq_len(mt_n_trait1_only)]
mt_trait2_only <- mt_active[mt_n_shared + mt_n_trait1_only + seq_len(mt_n_trait2_only)]

mt_state <- rep("00", nsnp)
mt_state[mt_trait1_only] <- "10"
mt_state[mt_trait2_only] <- "01"
mt_state[mt_shared] <- "11"

mt_beta1 <- numeric(nsnp)
mt_beta2 <- numeric(nsnp)

if (length(mt_trait1_only) > 0L) {
  freq_t1 <- colMeans(X[, mt_trait1_only, drop = FALSE]) / 2
  mt_beta1[mt_trait1_only] <- rnorm(length(mt_trait1_only), mean = 0, sd = sqrt((2 * freq_t1 * (1 - freq_t1))^S))
}
if (length(mt_trait2_only) > 0L) {
  freq_t2 <- colMeans(X[, mt_trait2_only, drop = FALSE]) / 2
  mt_beta2[mt_trait2_only] <- rnorm(length(mt_trait2_only), mean = 0, sd = sqrt((2 * freq_t2 * (1 - freq_t2))^S))
}
if (length(mt_shared) > 0L) {
  freq_shared <- colMeans(X[, mt_shared, drop = FALSE]) / 2
  shared_sd <- sqrt((2 * freq_shared * (1 - freq_shared))^S)
  z_shared_1 <- rnorm(length(mt_shared))
  z_shared_2 <- rnorm(length(mt_shared))
  mt_beta1[mt_shared] <- shared_sd * z_shared_1
  mt_beta2[mt_shared] <- shared_sd * (mt_shared_rho * z_shared_1 + sqrt(1 - mt_shared_rho^2) * z_shared_2)
}

g1 <- as.vector(X %*% mt_beta1)
g2 <- as.vector(X %*% mt_beta2)
g_mt <- cbind(g1, g2)
vg <- apply(g_mt, 2, var)
resid_var <- vg * (1 - h2) / h2
resid_cov <- mt_residual_rho * sqrt(resid_var[1] * resid_var[2])
resid_cov_matrix <- matrix(c(resid_var[1], resid_cov, resid_cov, resid_var[2]), 2, 2)
z_resid <- matrix(rnorm(nrow(X) * 2L), ncol = 2L)
e_mt <- z_resid %*% chol(resid_cov_matrix)
y_mt <- g_mt + e_mt

active_signal <- rnorm(nsnp, mean = ifelse(mt_state == "00", -0.3, 1.0), sd = 0.7)
pleiotropy_signal <- rnorm(nsnp, mean = ifelse(mt_state == "11", 1.0, ifelse(mt_state == "00", -0.2, 0.1)), sd = 0.7)
direction_signal <- rnorm(nsnp, mean = ifelse(mt_state == "10", 1.0, ifelse(mt_state == "01", -1.0, 0.0)), sd = 0.7)
random_signal <- rnorm(nsnp, mean = 0, sd = 1)

write.csv(
  data.frame(ID = ids, y1 = y_mt[1:n_gwas, 1], y2 = y_mt[1:n_gwas, 2]),
  file.path(dataset_dir, "phenotypes_mt.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    marker_id = marker_ids,
    active_signal = active_signal,
    pleiotropy_signal = pleiotropy_signal,
    direction_signal = direction_signal,
    random_signal = random_signal
  ),
  file.path(dataset_dir, "annotations_mt.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    marker_id = marker_ids,
    state = mt_state,
    is_active_y1 = mt_state %in% c("10", "11"),
    is_active_y2 = mt_state %in% c("01", "11"),
    is_shared = mt_state == "11",
    true_effect_y1 = mt_beta1,
    true_effect_y2 = mt_beta2
  ),
  file.path(dataset_dir, "truth_mt.csv"),
  row.names = FALSE
)

cat("Wrote simulated_annotations dataset to", dataset_dir, "\n")

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

cat("Wrote simulated_annotations dataset to", dataset_dir, "\n")

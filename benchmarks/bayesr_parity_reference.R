args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript benchmarks/bayesr_parity_reference.R <outdir> <fixed_pi|estimate_pi>")
}

outdir <- args[[1]]
mode <- args[[2]]
if (!(mode %in% c("fixed_pi", "estimate_pi"))) {
  stop("Mode must be fixed_pi or estimate_pi")
}

parse_config <- function(path) {
  cfg <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  setNames(cfg$value, cfg$key)
}

parse_numeric_vector <- function(value) {
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
}

rdirichlet1 <- function(alpha) {
  draws <- rgamma(length(alpha), shape = alpha, rate = 1)
  draws / sum(draws)
}

bayesr_reference <- function(X, y, niter, burnin, gamma, start_pi, start_sigma_sq, start_vare, estimate_pi) {
  n <- nrow(X)
  m <- ncol(X)
  ndist <- length(start_pi)
  if (niter <= burnin) {
    stop("niter must be greater than burnin")
  }

  pi <- start_pi
  sigma_sq <- start_sigma_sq
  vare <- start_vare
  nub <- 4
  nue <- 4
  scaleb <- (nub - 2) / nub * sigma_sq
  scalee <- (nue - 2) / nue * vare
  beta <- rep(0, m)
  mu <- mean(y)
  xpx <- colSums(X * X)
  ycorr <- y - mu

  kept_n <- niter - burnin
  beta_trace <- matrix(0, nrow = kept_n, ncol = m)
  pi_trace <- matrix(0, nrow = kept_n, ncol = ndist)
  sigma_trace <- numeric(kept_n)
  vare_trace <- numeric(kept_n)

  for (iter in seq_len(niter)) {
    ycorr <- ycorr + mu
    rhs_mu <- sum(ycorr)
    inv_lhs_mu <- 1 / n
    mu_hat <- inv_lhs_mu * rhs_mu
    mu <- rnorm(1, mu_hat, sqrt(inv_lhs_mu * vare))
    ycorr <- ycorr - mu

    log_pi <- log(pi)
    inv_sigma_sq <- rep(Inf, ndist)
    log_sigma_sq <- rep(Inf, ndist)
    nonzero_classes <- gamma > 0
    inv_sigma_sq[nonzero_classes] <- 1 / (gamma[nonzero_classes] * sigma_sq)
    log_sigma_sq[nonzero_classes] <- log(gamma[nonzero_classes] * sigma_sq)

    nsnp_dist <- rep(0, ndist)
    ssq <- 0

    for (j in seq_len(m)) {
      old_sample <- beta[j]
      rhs <- sum(X[, j] * ycorr) + xpx[j] * old_sample
      rhs <- rhs / vare
      inv_lhs <- 1 / (xpx[j] / vare + inv_sigma_sq)
      beta_hat <- inv_lhs * rhs

      log_delta <- 0.5 * (log(inv_lhs) - log_sigma_sq + beta_hat * rhs) + log_pi
      log_delta[1] <- log_pi[1]
      max_log <- max(log_delta)
      probs <- exp(log_delta - max_log)
      probs <- probs / sum(probs)

      delta <- sample.int(ndist, 1, prob = probs)
      nsnp_dist[delta] <- nsnp_dist[delta] + 1

      if (delta > 1) {
        beta[j] <- rnorm(1, beta_hat[delta], sqrt(inv_lhs[delta]))
        ycorr <- ycorr + X[, j] * (old_sample - beta[j])
        ssq <- ssq + beta[j]^2 / gamma[delta]
      } else {
        if (old_sample != 0) {
          ycorr <- ycorr + X[, j] * old_sample
        }
        beta[j] <- 0
      }
    }

    if (estimate_pi) {
      pi <- rdirichlet1(nsnp_dist + 1)
    }

    sigma_sq <- (ssq + nub * scaleb) / rchisq(1, sum(nsnp_dist[-1]) + nub)
    vare <- (sum(ycorr^2) + nue * scalee) / rchisq(1, n + nue)

    if (iter > burnin) {
      kept_idx <- iter - burnin
      beta_trace[kept_idx, ] <- beta
      pi_trace[kept_idx, ] <- pi
      sigma_trace[kept_idx] <- sigma_sq
      vare_trace[kept_idx] <- vare
    }
  }

  list(
    pi_mean = colMeans(pi_trace),
    sigma_sq_mean = mean(sigma_trace),
    vare_mean = mean(vare_trace),
    marker_estimate = colMeans(beta_trace),
    model_frequency = colMeans(beta_trace != 0)
  )
}

cfg <- parse_config(file.path(outdir, "data", "config.csv"))
seed <- as.integer(cfg[["seed"]])
niter <- as.integer(cfg[["chain_length"]])
burnin <- as.integer(cfg[["burnin"]])
gamma <- parse_numeric_vector(cfg[["gamma"]])
start_pi <- parse_numeric_vector(cfg[["start_pi"]])
start_sigma_sq <- as.numeric(cfg[["start_sigma_sq"]])
start_vare <- as.numeric(cfg[["start_vare"]])
estimate_pi <- identical(mode, "estimate_pi")

set.seed(seed)

geno_df <- read.csv(file.path(outdir, "data", "genotypes.csv"), stringsAsFactors = FALSE, check.names = FALSE)
marker_ids <- colnames(geno_df)[-1]
X <- as.matrix(geno_df[, -1, drop = FALSE])
storage.mode(X) <- "numeric"

pheno_df <- read.csv(file.path(outdir, "data", "phenotypes.csv"), stringsAsFactors = FALSE, check.names = FALSE)
y <- as.numeric(pheno_df[["y1"]])

fit <- bayesr_reference(X, y, niter, burnin, gamma, start_pi, start_sigma_sq, start_vare, estimate_pi)

summary_dir <- file.path(outdir, if (estimate_pi) "ref_estimate_pi" else "ref_fixed_pi")
if (dir.exists(summary_dir)) {
  unlink(summary_dir, recursive = TRUE, force = TRUE)
}
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  data.frame(
    metric = c("sigmaSq", "residual_variance", "mean_nonzero_frequency"),
    value = c(fit$sigma_sq_mean, fit$vare_mean, mean(fit$model_frequency))
  ),
  file.path(summary_dir, "scalar_metrics.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(
    class = paste0("class", seq_along(fit$pi_mean)),
    estimate = fit$pi_mean
  ),
  file.path(summary_dir, "pi.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(
    marker_id = marker_ids,
    estimate = fit$marker_estimate,
    model_frequency = fit$model_frequency
  ),
  file.path(summary_dir, "marker_effects.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("WROTE", summary_dir, "\n")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript benchmarks/bayesr_parity_replay_reference.R <outdir>")
}

outdir <- args[[1]]
datadir <- file.path(outdir, "data")

parse_config <- function(path) {
  cfg <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  setNames(cfg$value, cfg$key)
}

parse_numeric_vector <- function(value) {
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
}

parse_initial_scalars <- function(path) {
  scalars <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  setNames(as.numeric(scalars$value), scalars$key)
}

read_replay_draws <- function(path) {
  draws <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(colnames(draws), c("kind", "index", "value"))) {
    stop("replay draw file must contain kind,index,value columns")
  }
  draws
}

draw_value <- function(draws, kind, index = 1L) {
  rows <- draws[draws$kind == kind & draws$index == index, , drop = FALSE]
  if (nrow(rows) != 1) {
    stop(sprintf("Expected exactly one draw for %s[%d]", kind, index))
  }
  as.numeric(rows$value[[1]])
}

bayesr_class_probabilities <- function(rhs, xpx, vare, sigma_sq, pi, gamma) {
  nclasses <- length(gamma)
  log_probs <- rep(0, nclasses)
  log_probs[1] <- log(pi[1])
  inv_var_res <- 1 / vare
  for (k in 2:nclasses) {
    var_effect <- gamma[k] * sigma_sq
    inv_var_effect <- 1 / var_effect
    lhs <- xpx * inv_var_res + inv_var_effect
    inv_lhs <- 1 / lhs
    beta_hat <- inv_lhs * rhs
    log_probs[k] <- 0.5 * (log(inv_lhs) - log(var_effect) + beta_hat * rhs) + log(pi[k])
  }
  max_log <- max(log_probs)
  probs <- exp(log_probs - max_log)
  probs / sum(probs)
}

choose_class <- function(probs, u_class) {
  cumprob <- cumsum(probs)
  idx <- which(u_class <= cumprob)[1]
  if (is.na(idx)) {
    length(probs)
  } else {
    idx
  }
}

cfg <- parse_config(file.path(datadir, "config.csv"))
gamma <- parse_numeric_vector(cfg[["gamma"]])
pi <- parse_numeric_vector(cfg[["start_pi"]])
initial_scalars <- parse_initial_scalars(file.path(datadir, "initial_scalars.csv"))
draws <- read_replay_draws(file.path(datadir, "replay_draws_iteration1.csv"))

geno_df <- read.csv(file.path(datadir, "genotypes.csv"), stringsAsFactors = FALSE, check.names = FALSE)
marker_ids <- colnames(geno_df)[-1]
X <- as.matrix(geno_df[, -1, drop = FALSE])
storage.mode(X) <- "numeric"
xpx <- colSums(X * X)

pheno_df <- read.csv(file.path(datadir, "phenotypes.csv"), stringsAsFactors = FALSE, check.names = FALSE)
y <- as.numeric(pheno_df[["y1"]])
n <- length(y)

initial_state_df <- read.csv(file.path(datadir, "initial_state.csv"), stringsAsFactors = FALSE, check.names = FALSE)
beta <- as.numeric(initial_state_df$beta0)
delta <- as.integer(initial_state_df$delta0)
mu_old <- as.numeric(initial_scalars[["mu0"]])
sigma_sq_old <- as.numeric(initial_scalars[["sigmaSq0"]])
vare_old <- as.numeric(initial_scalars[["vare0"]])

ycorr <- as.numeric(y - mu_old - X %*% beta)
ycorr <- ycorr + mu_old
rhs_mu <- sum(ycorr)
inv_lhs_mu <- 1 / n
mu_hat <- inv_lhs_mu * rhs_mu
z_mu <- draw_value(draws, "mu_normal", 1L)
mu_new <- mu_hat + z_mu * sqrt(inv_lhs_mu * vare_old)
ycorr <- ycorr - mu_new

marker_rows <- vector("list", length(marker_ids))
inv_var_res <- 1 / vare_old

for (j in seq_along(marker_ids)) {
  old_alpha <- beta[j]
  rhs <- (sum(X[, j] * ycorr) + xpx[j] * old_alpha) * inv_var_res
  probs <- bayesr_class_probabilities(rhs, xpx[j], vare_old, sigma_sq_old, pi, gamma)
  u_class <- draw_value(draws, "marker_class_uniform", j)
  chosen_class <- choose_class(probs, u_class)
  delta[j] <- chosen_class

  beta_hat_chosen <- 0
  inv_lhs_chosen <- 0
  z_beta <- draw_value(draws, "marker_beta_normal", j)

  if (chosen_class == 1L) {
    if (old_alpha != 0) {
      ycorr <- ycorr + X[, j] * old_alpha
    }
    beta[j] <- 0
  } else {
    var_effect <- gamma[chosen_class] * sigma_sq_old
    inv_var_effect <- 1 / var_effect
    lhs <- xpx[j] * inv_var_res + inv_var_effect
    inv_lhs_chosen <- 1 / lhs
    beta_hat_chosen <- inv_lhs_chosen * rhs
    beta[j] <- beta_hat_chosen + z_beta * sqrt(inv_lhs_chosen)
    ycorr <- ycorr + X[, j] * (old_alpha - beta[j])
  }

  marker_rows[[j]] <- data.frame(
    marker_id = marker_ids[[j]],
    rhs = rhs,
    old_alpha = old_alpha,
    p_class1 = probs[[1]],
    p_class2 = probs[[2]],
    p_class3 = probs[[3]],
    p_class4 = probs[[4]],
    u_class = u_class,
    chosen_class = chosen_class,
    beta_hat_chosen = beta_hat_chosen,
    inv_lhs_chosen = inv_lhs_chosen,
    z_beta = z_beta,
    new_alpha = beta[j],
    ycorr_norm_after = sqrt(sum(ycorr^2))
  )
}

nonzero <- which(delta > 1L)
ssq <- if (length(nonzero) == 0) 0 else sum(beta[nonzero]^2 / gamma[delta[nonzero]])
nnz <- length(nonzero)
nub <- 4
nue <- 4
scaleb <- (nub - 2) / nub * sigma_sq_old
scalee <- (nue - 2) / nue * vare_old
chisq_sigma <- draw_value(draws, "sigma_chisq", 1L)
chisq_vare <- draw_value(draws, "vare_chisq", 1L)
sigma_sq_new <- (ssq + nub * scaleb) / chisq_sigma
vare_new <- (sum(ycorr^2) + nue * scalee) / chisq_vare

scalar_df <- data.frame(
  field = c(
    "mu_old",
    "mu_hat",
    "z_mu",
    "mu_new",
    "sigmaSq_old",
    "ssq",
    "nnz",
    "chisq_sigma",
    "sigmaSq_new",
    "vare_old",
    "chisq_vare",
    "vare_new"
  ),
  value = c(
    mu_old,
    mu_hat,
    z_mu,
    mu_new,
    sigma_sq_old,
    ssq,
    nnz,
    chisq_sigma,
    sigma_sq_new,
    vare_old,
    chisq_vare,
    vare_new
  )
)

summary_dir <- file.path(outdir, "ref_fixed_pi")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  do.call(rbind, marker_rows),
  file.path(summary_dir, "replay_marker_iteration1.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  scalar_df,
  file.path(summary_dir, "replay_scalar_iteration1.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("WROTE", summary_dir, "\n")

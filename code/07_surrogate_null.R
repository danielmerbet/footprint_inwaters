# 07_surrogate_null.R  -- REVISION S-D: selection-aware field-significance null
# Referee concern: BH-FDR treats pixels as independent and does not, by itself,
# establish FIELD significance against spatial+temporal autocorrelation; and the
# headline "median F among significant cells" is conditioned on significance.
# Here we build a phase-randomization (Ebisuzaki-type) null that PRESERVES each
# index's power spectrum (hence its autocorrelation) but destroys its true phase
# relationship to water temperature, rerun the WHOLE detection (per-cell adj-R^2,
# effective-DOF F-test, BH-FDR) on each surrogate, and compare the observed
# %-teleconnected and median-F to the null distribution.
#
# Engine is the lag-0 vectorised form of 03_footprint.R (justified: ~98% of cells
# are lag 0). Outputs: supplementary CSV + figure.
suppressMessages({library(ncdf4); library(terra); library(ggplot2); library(patchwork)})
## --- locate repo root (reproducible: Rscript, or source() in RStudio, from any dir) ---
local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  o <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  d <- if (length(a)) dirname(a[1]) else if (!is.null(o)) dirname(o) else getwd()
  d <- normalizePath(d, mustWork = FALSE)
  while (!all(dir.exists(file.path(d, c("input", "code")))) && dirname(d) != d) d <- dirname(d)
  assign("base", d, envir = .GlobalEnv)
})
setwd(base)
set.seed(42)
NSURR <- 500                       # number of surrogate index fields

years <- 1901:2021; keep <- years >= 1950 & years <= 2021; yrs <- years[keep]; n <- length(yrs)
X_t   <- cbind(1, scale(yrs, scale = FALSE))
Hcomp <- diag(n) - X_t %*% solve(crossprod(X_t)) %*% t(X_t)   # linear-detrend projector
detr  <- function(M) t(Hcomp %*% t(M))
si <- read.csv("input/clim_index/seasonal_indices.csv"); si <- si[match(yrs, si$year), ]
IDX0 <- scale(Hcomp %*% cbind(si$enso_djf, si$nao_djfm, si$iod_son))  # detrended, z-scored

systems <- list(
  surftemp = "input/laketemp/gotm_20crv3-era5_obsclim_histsoc_default_surftemp_global_annual_1901_2021.nc",
  bottemp  = "input/laketemp/gotm_20crv3-era5_obsclim_histsoc_default_bottemp_global_annual_1901_2021.nc",
  rivtemp  = "input/rivertemp/watergap2-2e_20crv3-era5_obsclim_histsoc_default_rivtemp_global_annual_1901_2021.nc")
vname <- c(surftemp = "surftemp", bottemp = "bottemp", rivtemp = "rivtemp")

# Ebisuzaki phase randomization: same power spectrum, random phases -> same
# autocorrelation, no true relationship to y.
phase_rand <- function(x) {
  N <- length(x); X <- fft(x); half <- if (N %% 2 == 0) N/2 else (N + 1)/2
  k <- 2:half; phi <- runif(length(k), 0, 2 * pi)
  X[k] <- Mod(X[k]) * exp(1i * phi)
  X[N - k + 2] <- Conj(X[k])                     # enforce conjugate symmetry -> real output
  Re(fft(X, inverse = TRUE) / N)
}

# Vectorised detection over all cells for a given (detrended, z-scored) index matrix.
# Returns per-cell adjusted R^2 and effective-DOF BH q-value (the full 03 engine at lag 0).
detect <- function(Yd, IDX) {
  p <- ncol(IDX) + 1; X <- cbind(1, IDX)
  H <- X %*% solve(crossprod(X)) %*% t(X)        # n x n hat matrix (common to all cells)
  res    <- Yd - Yd %*% t(H)                     # np x n
  sstot  <- rowSums((Yd - rowMeans(Yd))^2)
  r2  <- 1 - rowSums(res^2) / sstot
  adj <- 1 - (1 - r2) * (n - 1) / (n - p)
  a <- res[, -1, drop = FALSE]; b <- res[, -n, drop = FALSE]     # lag-1 autocorr of residuals
  am <- a - rowMeans(a); bm <- b - rowMeans(b)
  r1 <- rowSums(am * bm) / sqrt(rowSums(am^2) * rowSums(bm^2))
  r1[is.na(r1)] <- 0; r1 <- pmin(pmax(r1, 0), 0.99)
  neff <- n * (1 - r1) / (1 + r1)
  Fst  <- (r2 / (p - 1)) / ((1 - r2) / (neff - p))
  pv <- rep(NA_real_, length(r2)); ok <- neff - p > 1 & Fst > 0 & sstot > 0 & is.finite(Fst)
  pv[ok] <- pf(Fst[ok], p - 1, neff[ok] - p, lower.tail = FALSE)
  list(adj = adj, q = p.adjust(pv, method = "BH"))
}

outdir <- "output/figures/supplementary/robustness/"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
summ <- list(); nulls <- list()
for (nm in names(systems)) {
  cat("###### surrogate null:", nm, "######\n")
  nc <- nc_open(systems[[nm]]); v <- ncvar_get(nc, vname[nm]); nc_close(nc); v <- v[, , keep]
  grid <- rast(lapply(seq_len(n), function(k){r<-rast(t(v[,,k]));ext(r)<-c(-180,180,-90,90);crs(r)<-"EPSG:4326";r}))
  tab <- as.data.frame(grid, xy = TRUE, na.rm = TRUE); Yd <- detr(as.matrix(tab[, -(1:2)]))

  obs <- detect(Yd, IDX0); npix <- length(obs$adj)
  sigobs <- which(obs$q < 0.05 & !is.na(obs$q)); Nobs <- length(sigobs)
  obs_pct <- 100 * Nobs / npix; obs_medF <- median(obs$adj[sigobs], na.rm = TRUE)
  # null: (a) extent = % passing FDR; (b) selection-aware strength = median F of the
  # top-Nobs cells drawn from the noise field (tests the winner's curse directly).
  null_pct <- numeric(NSURR); null_topF <- numeric(NSURR)
  for (s in seq_len(NSURR)) {
    IDXs <- scale(apply(IDX0, 2, phase_rand))     # randomize each index independently
    d <- detect(Yd, IDXs)
    null_pct[s]  <- 100 * sum(d$q < 0.05 & !is.na(d$q)) / npix
    null_topF[s] <- median(sort(d$adj, decreasing = TRUE)[seq_len(Nobs)], na.rm = TRUE)
  }
  p_pct <- (1 + sum(null_pct  >= obs_pct )) / (NSURR + 1)
  p_str <- (1 + sum(null_topF >= obs_medF)) / (NSURR + 1)
  cat(sprintf("  obs %%-teleconnected = %.2f | null median = %.2f [95%% = %.2f] | p = %.4f\n",
              obs_pct, median(null_pct), quantile(null_pct, .95), p_pct))
  cat(sprintf("  obs median F = %.3f | null top-%d median F = %.3f [95%% = %.3f] | p = %.4f\n",
              obs_medF, Nobs, median(null_topF), quantile(null_topF, .95), p_str))
  summ[[nm]] <- data.frame(system = nm, n_sig = Nobs, obs_pct = obs_pct, null_pct_med = median(null_pct),
    null_pct_p95 = quantile(null_pct, .95), p_pct = p_pct,
    obs_medF = obs_medF, null_topF_med = median(null_topF),
    null_topF_p95 = quantile(null_topF, .95), p_strength = p_str, nsurr = NSURR)
  nulls[[nm]] <- data.frame(system = nm, pct = null_pct, topF = null_topF)
}
S <- do.call(rbind, summ); rownames(S) <- NULL
write.csv(S, paste0(outdir, "surrogate_null.csv"), row.names = FALSE)

ND <- do.call(rbind, nulls)
lab <- c(surftemp="Lake surface", bottemp="Lake bottom", rivtemp="River")
ND$sys <- factor(lab[ND$system], levels = lab); S$sys <- factor(lab[S$system], levels = lab)
pA <- ggplot(ND, aes(pct)) + geom_histogram(bins = 40, fill = "grey75", colour = "white") +
  geom_vline(data = S, aes(xintercept = obs_pct), colour = "#C1272D", linewidth = 1) +
  facet_wrap(~sys, scales = "free") +
  labs(title = "a  Field significance: % teleconnected vs phase-randomized null",
       subtitle = "grey = null (spectrum-preserving surrogates); red = observed",
       x = "% of cells teleconnected", y = "surrogate count") +
  theme_minimal(base_size = 11)
pB <- ggplot(ND, aes(topF)) + geom_histogram(bins = 40, fill = "grey75", colour = "white") +
  geom_vline(data = S, aes(xintercept = obs_medF), colour = "#C1272D", linewidth = 1) +
  facet_wrap(~sys, scales = "free") +
  labs(title = "b  Selection-aware strength: observed median F vs null median F of the top-N cells",
       subtitle = "grey = median F of the N best cells drawn from noise; red = observed median F",
       x = "median F (adj. R²)", y = "surrogate count") +
  theme_minimal(base_size = 11)
ggsave(paste0(outdir, "surrogate_null.png"), pA / pB, width = 10, height = 6.4, dpi = 200, bg = "white")
cat("\nsaved -> ", outdir, "surrogate_null.{csv,png}\n", sep = "")

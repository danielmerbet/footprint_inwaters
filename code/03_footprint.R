# 03_footprint.R  -- oscillation footprint F for every inland-water compartment
# Fraction of INTERANNUAL water-temperature variance organised by ENSO, IOD and NAO,
# fully self-contained & detrended-consistent. One parameterised engine for all three
# compartments (lake surface, lake bottom, river).
#
# Per pixel, on the annual series 1950-2021:
#   1. detrend water temp AND indices (linear -> interannual band).
#   2. season-specific indices (ENSO=DJF, NAO=DJFM, IOD=SON), z-scored.
#   3. per-index detrended significance -> coherent patches (>= min_region) -> cluster lag.
#   4. footprint F = adjusted R^2 of  y' ~ (lagged detrended indices)   [ = R^2_O ].
#   5. effective-DOF significance + BH-FDR field significance.
# Outputs -> output/footprint/{compartment}/ (GeoTIFF).

suppressMessages({library(ncdf4); library(terra)})
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
indices    <- c("enso", "nao", "iod")
maxlag     <- 3
min_region <- 5
index_col  <- c(enso = "enso_djf", nao = "nao_djfm", iod = "iod_son")

systems <- list(
  surftemp = list(nc = "input/laketemp/gotm_20crv3-era5_obsclim_histsoc_default_surftemp_global_annual_1901_2021.nc",       var = "surftemp"),
  bottemp  = list(nc = "input/laketemp/gotm_20crv3-era5_obsclim_histsoc_default_bottemp_global_annual_1901_2021.nc",        var = "bottemp"),
  rivtemp  = list(nc = "input/rivertemp/watergap2-2e_20crv3-era5_obsclim_histsoc_default_rivtemp_global_annual_1901_2021.nc", var = "rivtemp"))

run_footprint <- function(name, cfg) {
  cat("###### footprint:", name, "######\n")
  nc <- nc_open(cfg$nc); v <- ncvar_get(nc, cfg$var); nc_close(nc)
  years <- 1901:2021; keep <- years >= 1950 & years <= 2021
  v <- v[, , keep]; yrs <- years[keep]; n <- length(yrs)
  grid <- rast(lapply(seq_len(n), function(k) {
    r <- rast(t(v[, , k])); ext(r) <- c(-180, 180, -90, 90); crs(r) <- "EPSG:4326"; r }))
  tab <- as.data.frame(grid, xy = TRUE, na.rm = TRUE)
  Y <- as.matrix(tab[, -(1:2)]); xy <- tab[, 1:2]
  cells <- cellFromXY(grid[[1]], as.matrix(xy)); np <- nrow(Y)
  cat("  pixels:", np, " years:", n, "\n")

  X_t <- cbind(1, scale(yrs, scale = FALSE))
  Hcomp <- diag(n) - X_t %*% solve(crossprod(X_t)) %*% t(X_t)
  Yd <- t(Hcomp %*% t(Y))
  si <- read.csv("input/clim_index/seasonal_indices.csv"); si <- si[match(yrs, si$year), ]
  IDX <- sapply(indices, function(i) si[[index_col[i]]])
  IDXd <- scale(Hcomp %*% IDX)

  pixel_corr_p <- function(ind) {
    yc <- Yd - rowMeans(Yd); ic <- ind - mean(ind)
    r <- as.numeric((yc %*% ic) / (sqrt(rowSums(yc^2)) * sqrt(sum(ic^2))))
    list(r = r, p = 2 * pt(-abs(r * sqrt((n - 2) / (1 - r^2))), df = n - 2)) }
  best_lag_cor <- function(series, ind) {
    cors <- sapply(0:maxlag, function(L) suppressWarnings(cor(series[(1 + L):n], ind[1:(n - L)], use = "complete.obs")))
    (0:maxlag)[which.max(abs(cors))] }
  LAGS <- matrix(0L, np, length(indices), dimnames = list(NULL, indices))
  SIG  <- matrix(FALSE, np, length(indices), dimnames = list(NULL, indices))
  for (j in seq_along(indices)) {
    cp <- pixel_corr_p(IDXd[, j]); sig <- cp$p < 0.05 & !is.na(cp$p); SIG[, j] <- sig
    mr <- grid[[1]]; values(mr) <- NA; mr[cells[sig]] <- 1
    reg <- patches(mr, directions = 8, zeroAsNA = TRUE)
    rid <- as.numeric(terra::extract(reg, as.matrix(xy))[, 1])
    for (g in unique(rid[!is.na(rid)])) {
      rows <- which(rid == g); if (length(rows) < min_region) next
      LAGS[rows, j] <- best_lag_cor(colMeans(Yd[rows, , drop = FALSE], na.rm = TRUE), IDXd[, j]) } }

  fit_adjR2 <- function(y, X) {
    p <- ncol(X); nn <- length(y); fit <- .lm.fit(X, y); res <- fit$residuals
    ss_res <- sum(res^2); ss_tot <- sum((y - mean(y))^2)
    if (ss_tot <= 0 || nn - p <= 1) return(c(NA, NA))
    r2 <- 1 - ss_res / ss_tot; adj <- 1 - (1 - r2) * (nn - 1) / (nn - p)
    r1 <- suppressWarnings(cor(res[-1], res[-nn])); r1 <- ifelse(is.na(r1), 0, max(min(r1, 0.99), 0))
    neff <- nn * (1 - r1) / (1 + r1); Fst <- (r2 / (p - 1)) / ((1 - r2) / (neff - p))
    pv <- if (neff - p > 1 && Fst > 0) pf(Fst, p - 1, neff - p, lower.tail = FALSE) else NA
    c(adj, pv) }
  F_adj <- rep(NA_real_, np); p_eff <- rep(NA_real_, np)
  adj_single <- matrix(NA_real_, np, length(indices), dimnames = list(NULL, indices))
  for (i in seq_len(np)) {
    L <- LAGS[i, ]; mL <- max(L); rows <- (1 + mL):n; y <- Yd[i, rows]
    Xi <- sapply(seq_along(indices), function(j) IDXd[rows - L[j], j])
    o <- fit_adjR2(y, cbind(1, Xi)); F_adj[i] <- o[1]; p_eff[i] <- o[2]
    for (j in seq_along(indices)) adj_single[i, j] <- fit_adjR2(y, cbind(1, IDXd[rows - L[j], j]))[1] }
  q_fdr <- p.adjust(p_eff, method = "BH")

  outdir <- paste0("output/footprint/", name, "/"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  wm <- function(vals, nm) { r <- grid[[1]]; values(r) <- NA; r[cells] <- vals; writeRaster(r, paste0(outdir, nm, ".tiff"), overwrite = TRUE) }
  wm(F_adj, "footprint_adjR2"); wm(p_eff, "pvalue_eff"); wm(q_fdr, "qvalue_fdr")
  for (j in seq_along(indices)) { wm(adj_single[, j], paste0(indices[j], "_adjR2"))
    wm(LAGS[, j], paste0(indices[j], "_lag")); wm(as.integer(SIG[, j]), paste0(indices[j], "_sig")) }
  sigfdr <- q_fdr < 0.05 & !is.na(q_fdr)
  cat(sprintf("  median F (FDR-sig) = %.3f | %% teleconnected = %.1f\n\n",
              median(F_adj[sigfdr], na.rm = TRUE), 100 * sum(sigfdr) / np))
}
for (nm in names(systems)) run_footprint(nm, systems[[nm]])
cat("done -> output/footprint/{surftemp,bottemp,rivtemp}/\n")

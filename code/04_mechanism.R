# 04_mechanism.R  -- how the oscillation signal reaches each teleconnected cell
# Two-block variation partitioning (Borcard/Legendre; a+b = F) splitting the footprint into:
#   direct   [a] = oscillation variance NOT shared with local climate (remote / thermal memory)
#   mediated [b] = shared oscillation & local climate (delivered through local weather)
#   local    [c] = local-climate variance not organised by the oscillations
# Rivers add PRECIPITATION to the local block (hydrological pathway lakes lack).
# Outputs -> output/footprint/{compartment}/mechanism/

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
years <- 1901:2021; keep <- years >= 1950 & years <= 2021; yrs <- years[keep]; n <- length(yrs)
X_t <- cbind(1, scale(yrs, scale = FALSE)); Hcomp <- diag(n) - X_t %*% solve(crossprod(X_t)) %*% t(X_t)
detmat <- function(M) t(Hcomp %*% t(M))
si <- read.csv("input/clim_index/seasonal_indices.csv"); si <- si[match(yrs, si$year), ]
IDXd <- scale(Hcomp %*% cbind(si$enso_djf, si$nao_djfm, si$iod_son)); E <- IDXd[,1]; N <- IDXd[,2]; I <- IDXd[,3]

systems <- list(
  surftemp = list(nc = "input/laketemp/gotm_20crv3-era5_obsclim_histsoc_default_surftemp_global_annual_1901_2021.nc",       var = "surftemp", locvars = c("tas","rsds","rlds","sfcwind")),
  bottemp  = list(nc = "input/laketemp/gotm_20crv3-era5_obsclim_histsoc_default_bottemp_global_annual_1901_2021.nc",        var = "bottemp",  locvars = c("tas","rsds","rlds","sfcwind")),
  rivtemp  = list(nc = "input/rivertemp/watergap2-2e_20crv3-era5_obsclim_histsoc_default_rivtemp_global_annual_1901_2021.nc", var = "rivtemp",  locvars = c("tas","rsds","rlds","sfcwind","pr")))

adjR2 <- function(y, X) { p <- ncol(X); nn <- length(y); f <- .lm.fit(X, y)
  ssr <- sum(f$residuals^2); sst <- sum((y - mean(y))^2)
  if (sst <= 0 || nn - p <= 1) return(NA); 1 - (1 - (1 - ssr / sst)) * (nn - 1) / (nn - p) }

run_mech <- function(name, cfg) {
  cat("###### mechanism:", name, "| local:", paste(cfg$locvars, collapse=","), "######\n")
  nc <- nc_open(cfg$nc); v <- ncvar_get(nc, cfg$var); nc_close(nc); v <- v[, , keep]
  grid <- rast(lapply(seq_len(n), function(k){r<-rast(t(v[,,k]));ext(r)<-c(-180,180,-90,90);crs(r)<-"EPSG:4326";r}))
  tab <- as.data.frame(grid, xy = TRUE, na.rm = TRUE); xy <- tab[, 1:2]
  cells <- cellFromXY(grid[[1]], as.matrix(xy)); Yd <- detmat(as.matrix(tab[, -(1:2)]))
  d <- paste0("output/footprint/", name, "/"); gv <- function(f) rast(paste0(d, f))[cells][, 1]
  q <- gv("qvalue_fdr.tiff"); lE <- gv("enso_lag.tiff"); lN <- gv("nao_lag.tiff"); lI <- gv("iod_lag.tiff")
  loadmet <- function(var) { nc <- nc_open(paste0("input/clim_data/20crv3-era5_obsclim_", var, "_global_yearly_1901_2021.nc"))
    m <- ncvar_get(nc, var); nc_close(nc); m <- m[, , keep]
    st <- rast(lapply(seq_len(n), function(k){r<-rast(t(m[,,k]));ext(r)<-c(-180,180,-90,90);crs(r)<-"EPSG:4326";r})); detmat(st[cells]) }
  L <- lapply(cfg$locvars, loadmet); names(L) <- cfg$locvars

  np <- nrow(Yd); direct<-mediated<-local<-medfrac<-rep(NA_real_,np); dommed<-rep(NA_integer_,np)
  for (i in which(q < 0.05 & !is.na(q))) {
    Lg <- c(lE[i], lN[i], lI[i]); mL <- max(Lg); rows <- (1 + mL):n; y <- Yd[i, rows]
    O <- cbind(E[rows-lE[i]], N[rows-lN[i]], I[rows-lI[i]]); Lb <- sapply(cfg$locvars, function(vn) L[[vn]][i, rows])
    R2o <- adjR2(y, cbind(1, O)); R2l <- adjR2(y, cbind(1, Lb)); R2f <- adjR2(y, cbind(1, O, Lb))
    if (any(is.na(c(R2o, R2l, R2f)))) next
    a <- R2f - R2l; b <- R2o - a; direct[i]<-a; mediated[i]<-b; local[i]<-R2f-R2o
    medfrac[i] <- if (R2o > 0.01) max(0, min(1, b / R2o)) else NA
    dommed[i] <- which.max(sapply(cfg$locvars, function(vn) cor(L[[vn]][i, rows], fitted(lm(y ~ O)))^2)) }

  outdir <- paste0(d, "mechanism/"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  wr <- function(vals, nm){r<-grid[[1]];values(r)<-NA;r[cells]<-vals;writeRaster(r,paste0(outdir,nm,".tiff"),overwrite=TRUE)}
  wr(direct,"direct_adjR2"); wr(mediated,"mediated_adjR2"); wr(local,"local_adjR2"); wr(medfrac,"mediated_frac"); wr(dommed,"dominant_mediator")
  ok <- !is.na(medfrac)
  cat(sprintf("  cells=%d | median mediated frac=%.2f | dominant mediator: %s\n\n", sum(ok), median(medfrac[ok],na.rm=TRUE),
      paste(sprintf("%s %.0f%%", cfg$locvars, 100*prop.table(table(factor(dommed[ok],seq_along(cfg$locvars))))), collapse=" ")))
}
for (nm in names(systems)) run_mech(nm, systems[[nm]])
cat("done -> output/footprint/{...}/mechanism/\n")

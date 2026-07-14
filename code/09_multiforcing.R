# 09_multiforcing.R  -- REVISION S-F: is the footprint a single-forcing artifact?
# Referee concern: lake temperatures are one model (GOTM) driven by one reanalysis.
# We cannot swap the lake model here, but the ISIMIP GOTM runs are available under
# THREE meteorological forcings (20CRv3-ERA5, 20CRv3-W5E5, GSWP3-W5E5). If the
# oscillation footprint is physical it should replicate across forcings; if it is an
# artifact of one reanalysis it should not. We recompute F (lag-0 vectorised engine,
# common period 1950-2019) for lake surface and bottom under each forcing and compare
# the F fields (pixelwise correlation), the median F, and the teleconnected-cell overlap.
# Rivers have a single forcing (WaterGAP2 triver) and are excluded from this test.
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
CY <- 1950:2019                       # common period across forcings

si <- read.csv("input/clim_index/seasonal_indices.csv")

# lag-0 vectorised footprint -> full-grid rasters of adj-R^2 and BH q-value
footprint_field <- function(ncfile, var, fileyears) {
  keep <- fileyears %in% CY; yrs <- fileyears[keep]; n <- length(yrs)
  nc <- nc_open(ncfile); v <- ncvar_get(nc, var); nc_close(nc); v <- v[, , keep]
  grid <- rast(lapply(seq_len(n), function(k){r<-rast(t(v[,,k]));ext(r)<-c(-180,180,-90,90);crs(r)<-"EPSG:4326";r}))
  tab <- as.data.frame(grid, xy = TRUE, na.rm = TRUE)
  cells <- cellFromXY(grid[[1]], as.matrix(tab[, 1:2]))
  X_t <- cbind(1, scale(yrs, scale = FALSE)); H0 <- diag(n) - X_t %*% solve(crossprod(X_t)) %*% t(X_t)
  Yd <- t(H0 %*% t(as.matrix(tab[, -(1:2)])))
  s  <- si[match(yrs, si$year), ]
  IDX <- scale(H0 %*% cbind(s$enso_djf, s$nao_djfm, s$iod_son))
  p <- 4; X <- cbind(1, IDX); Hh <- X %*% solve(crossprod(X)) %*% t(X)
  res <- Yd - Yd %*% t(Hh)
  sstot <- rowSums((Yd - rowMeans(Yd))^2); r2 <- 1 - rowSums(res^2)/sstot
  adj <- 1 - (1 - r2) * (n - 1)/(n - p)
  a <- res[,-1,drop=FALSE]; b <- res[,-n,drop=FALSE]; am<-a-rowMeans(a); bm<-b-rowMeans(b)
  r1 <- rowSums(am*bm)/sqrt(rowSums(am^2)*rowSums(bm^2)); r1[is.na(r1)]<-0; r1<-pmin(pmax(r1,0),0.99)
  neff <- n*(1-r1)/(1+r1); Fst <- (r2/(p-1))/((1-r2)/(neff-p))
  pv <- rep(NA_real_, length(r2)); ok <- neff-p>1 & Fst>0 & sstot>0 & is.finite(Fst)
  pv[ok] <- pf(Fst[ok], p-1, neff[ok]-p, lower.tail=FALSE); q <- p.adjust(pv, "BH")
  radj <- grid[[1]]; values(radj) <- NA; radj[cells] <- adj
  rq   <- grid[[1]]; values(rq)   <- NA; rq[cells]   <- q
  list(F = radj, q = rq)
}

forcings <- list(
  `20CRv3-ERA5` = list(sfx="20crv3-era5",  yr=1901:2021),
  `20CRv3-W5E5` = list(sfx="20crv3-w5e5",  yr=1901:2019),
  `GSWP3-W5E5`  = list(sfx="gswp3-w5e5",   yr=1901:2019))
compartments <- c(surftemp = "surftemp", bottemp = "bottemp")
lab <- c(surftemp = "Lake surface", bottemp = "Lake bottom")

outdir <- "output/figures/supplementary/robustness/"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
rows <- list(); plots <- list()
for (cp in names(compartments)) {
  cat("\n########## multi-forcing:", lab[cp], "(", paste(range(CY),collapse="-"), ") ##########\n")
  fields <- lapply(forcings, function(f)
    footprint_field(sprintf("input/laketemp/gotm_%s_obsclim_histsoc_default_%s_global_annual_%d_%d.nc",
                            f$sfx, compartments[cp], min(f$yr), max(f$yr)), compartments[cp], f$yr))
  base_f <- fields[["20CRv3-ERA5"]]
  for (fn in names(forcings)) {
    Ff <- fields[[fn]]$F; qf <- fields[[fn]]$q
    sig <- values(qf)[,1] < 0.05 & !is.na(values(qf)[,1])
    medF <- median(values(Ff)[,1][sig], na.rm = TRUE); pct <- 100*sum(sig)/sum(!is.na(values(Ff)[,1]))
    # agreement vs ERA5 baseline
    both <- !is.na(values(base_f$F)[,1]) & !is.na(values(Ff)[,1])
    rF <- cor(values(base_f$F)[,1][both], values(Ff)[,1][both])
    bsig <- values(base_f$q)[,1] < 0.05 & !is.na(values(base_f$q)[,1])
    jac <- sum(sig & bsig)/sum(sig | bsig)
    cat(sprintf("  %-12s median F=%.3f | %%tele=%.2f | pixel r(F) vs ERA5=%.3f | Jaccard(sig) vs ERA5=%.2f\n",
                fn, medF, pct, rF, jac))
    rows[[length(rows)+1]] <- data.frame(compartment=cp, forcing=fn, medF=medF, pct=pct, rF=rF, jaccard=jac)
  }
  # scatter of ERA5 vs each alternative (sig-in-either)
  for (fn in c("20CRv3-W5E5","GSWP3-W5E5")) {
    both <- !is.na(values(base_f$F)[,1]) & !is.na(values(fields[[fn]]$F)[,1])
    dd <- data.frame(era5=values(base_f$F)[,1][both], alt=values(fields[[fn]]$F)[,1][both])
    plots[[paste(cp,fn)]] <- ggplot(dd, aes(era5, alt)) +
      geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey60") +
      geom_point(alpha=.05, size=.5, colour="#3b528b") +
      coord_equal(xlim=c(0,.6), ylim=c(0,.6)) +
      labs(title=sprintf("%s: ERA5 vs %s", lab[cp], fn), x="F (20CRv3-ERA5)", y=paste0("F (",fn,")")) +
      theme_minimal(base_size=10)
  }
}
M <- do.call(rbind, rows); write.csv(M, paste0(outdir, "multiforcing.csv"), row.names = FALSE)
ggsave(paste0(outdir, "multiforcing.png"), wrap_plots(plots, ncol = 2), width = 9, height = 8, dpi = 200, bg = "white")
cat("\nsaved -> ", outdir, "multiforcing.{csv,png}\n", sep = "")

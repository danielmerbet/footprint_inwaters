# 08_depth_memory.R  -- REVISION S-E: does the direct/memory pathway scale with lake depth?
# The central "buffered deep lake" claim predicts that the DIRECT (non-locally-mediated)
# fraction of the footprint should grow with lake depth. We test this directly using the
# previously unused bathymetry (input/lake_bathymetry/max_lakedepth.nc): for teleconnected
# lake cells, regress direct fraction (= 1 - mediated fraction) on maximum lake depth,
# separately for surface and bottom. A positive depth relationship at the bottom (and its
# absence at the surface) is independent, morphometric support for the memory interpretation.
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

# --- load bathymetry, align to the footprint grid ---
dpath <- "input/lake_bathymetry/max_lakedepth.nc"
dnc <- nc_open(dpath)
dvar <- setdiff(names(dnc$var), c("lon","lat","time","crs","lon_bnds","lat_bnds"))
if (length(dvar) != 1) dvar <- dvar[which.max(sapply(dvar, function(vn) prod(dnc$var[[vn]]$size)))]
dm <- ncvar_get(dnc, dvar); nc_close(dnc)
depth <- rast(t(dm)); ext(depth) <- c(-180,180,-90,90); crs(depth) <- "EPSG:4326"
fpref <- rast("output/footprint/bottemp/footprint_adjR2.tiff")
if (!isTRUE(compareGeom(depth, fpref, stopOnError = FALSE)))
  depth <- resample(depth, fpref, method = "bilinear")
dv <- values(depth)[,1]
cat(sprintf("depth var='%s' | range = %.1f - %.1f m | median = %.1f m\n",
            dvar, min(dv,na.rm=TRUE), max(dv,na.rm=TRUE), median(dv,na.rm=TRUE)))

grab <- function(compartment) {
  d  <- paste0("output/footprint/", compartment, "/")
  q  <- rast(paste0(d, "qvalue_fdr.tiff"))
  mf <- rast(paste0(d, "mechanism/mediated_frac.tiff"))
  sig <- q < 0.05
  df <- data.frame(depth = values(depth)[,1], medfrac = values(mf)[,1], sig = values(sig)[,1])
  df <- df[which(df$sig == 1 & !is.na(df$medfrac) & !is.na(df$depth) & df$depth > 0), ]
  df$direct <- pmax(0, pmin(1, 1 - df$medfrac))
  df$compartment <- compartment; df
}
S <- grab("surftemp"); B <- grab("bottemp")
lab <- c(surftemp = "Lake surface", bottemp = "Lake bottom")

report <- function(df, nm) {
  rho <- suppressWarnings(cor(log10(df$depth), df$direct, method = "spearman"))
  fit <- lm(direct ~ log10(depth), data = df)
  cat(sprintf("%-12s n=%5d | Spearman(direct, log10 depth) = %+.3f | slope/decade = %+.3f (p=%.2g)\n",
              nm, nrow(df), rho, coef(fit)[2], summary(fit)$coefficients[2,4]))
  # binned medians for the paper text
  br <- c(0,3,10,30,100,Inf); df$bin <- cut(df$depth, br, labels=c("<3","3-10","10-30","30-100",">100"))
  bm <- aggregate(direct ~ bin, df, median)
  cat("   median direct fraction by depth (m):",
      paste(sprintf("%s=%.3f", bm$bin, bm$direct), collapse="  "), "\n")
  invisible(rho)
}
cat("\n--- direct (memory) fraction vs maximum lake depth ---\n")
report(S, lab["surftemp"]); report(B, lab["bottemp"])

D <- rbind(S, B); D$sys <- factor(lab[D$compartment], levels = lab)
outdir <- "output/figures/supplementary/robustness/"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
p <- ggplot(D, aes(depth, direct)) +
  geom_point(alpha = .12, size = .6, colour = "#2E86AB") +
  geom_smooth(method = "loess", se = TRUE, colour = "#C1272D", fill = "#C1272D22", linewidth = .9) +
  scale_x_log10() + coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~sys) +
  labs(title = "Direct (memory) fraction of the oscillation footprint vs lake depth",
       subtitle = "teleconnected cells; direct = 1 - mediated fraction",
       x = "maximum lake depth (m, log scale)", y = "direct fraction of F") +
  theme_minimal(base_size = 11)
ggsave(paste0(outdir, "depth_memory.png"), p, width = 9, height = 4.2, dpi = 200, bg = "white")
write.csv(D[, c("compartment","depth","direct")], paste0(outdir, "depth_memory.csv"), row.names = FALSE)
cat("\nsaved -> ", outdir, "depth_memory.{csv,png}\n", sep = "")

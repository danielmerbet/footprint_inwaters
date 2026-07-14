# export_interactive_data.R  (RIVERS)
# Compact JSON for the interactive river-sentinel explorer.
#  - shared detrended z-scored season indices E/N/I (72 yr)
#  - per FDR-significant river cell: detrended obs series + model coefs/lags + per-index R2
#  - all river cells: lon/lat + F + dominant index, for the grey basemap & click info
suppressMessages({library(ncdf4); library(terra); library(jsonlite)})
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
X_t   <- cbind(1, scale(yrs, scale = FALSE))
Hcomp <- diag(n) - X_t %*% solve(crossprod(X_t)) %*% t(X_t)

si  <- read.csv("input/clim_index/seasonal_indices.csv"); si <- si[match(yrs, si$year), ]
IDX  <- cbind(si$enso_djf, si$nao_djfm, si$iod_son)
IDXd <- scale(Hcomp %*% IDX)                      # cols: ENSO, NAO, IOD
E <- IDXd[,1]; N <- IDXd[,2]; I <- IDXd[,3]

variable <- "rivtemp"
nc <- nc_open(paste0("input/rivertemp/watergap2-2e_20crv3-era5_obsclim_histsoc_default_",
                     variable, "_global_annual_1901_2021.nc"))
v <- ncvar_get(nc, variable); nc_close(nc); v <- v[, , keep]
lake <- rast(lapply(seq_len(n), function(k){r<-rast(t(v[,,k]));ext(r)<-c(-180,180,-90,90);crs(r)<-"EPSG:4326";r}))
tab  <- as.data.frame(lake, xy = TRUE, na.rm = TRUE)
Y    <- as.matrix(tab[,-(1:2)]); xy <- tab[,1:2]
Yd   <- t(Hcomp %*% t(Y))
d    <- paste0("output/footprint/", variable, "/")
cells<- cellFromXY(lake[[1]], as.matrix(xy))
getv <- function(f) rast(paste0(d,f))[cells][,1]
V <- list(xy=xy, Yd=Yd,
     F=getv("footprint_adjR2.tiff"), q=getv("qvalue_fdr.tiff"),
     lE=getv("enso_lag.tiff"), lN=getv("nao_lag.tiff"), lI=getv("iod_lag.tiff"),
     r2e=getv("enso_adjR2.tiff"), r2n=getv("nao_adjR2.tiff"), r2i=getv("iod_adjR2.tiff"))

sig <- which(V$q < 0.05 & !is.na(V$q))
recs <- lapply(sig, function(i) {
  L <- c(V$lE[i], V$lN[i], V$lI[i]); mL <- max(L); rows <- (1+mL):n
  Xi <- cbind(1, E[rows-V$lE[i]], N[rows-V$lN[i]], I[rows-V$lI[i]])
  b  <- coef(.lm.fit(Xi, V$Yd[i, rows]))
  list(lon=round(V$xy$x[i],3), lat=round(V$xy$y[i],3),
       F=round(V$F[i],3), dom=which.max(c(V$r2e[i],V$r2i[i],V$r2n[i])),  # 1=ENSO 2=IOD 3=NAO
       r2=round(c(V$r2e[i],V$r2i[i],V$r2n[i]),2),
       a=round(b[1],3), cE=round(b[2],3), cN=round(b[3],3), cI=round(b[4],3),
       lag=c(V$lE[i],V$lN[i],V$lI[i]),
       obs=round(V$Yd[i,],2))
})

dom_all <- apply(cbind(V$r2e,V$r2i,V$r2n),1,which.max); dom_all[V$q>=0.05|is.na(V$q)] <- 0
ctx <- list(lon=round(V$xy$x,3), lat=round(V$xy$y,3), F=round(V$F,2), dom=dom_all)

out <- list(
  years = yrs,
  idx   = list(E=round(E,3), N=round(N,3), I=round(I,3)),
  dom_names = c("ENSO","IOD","NAO"),
  sig   = recs,
  ctx   = ctx
)
js <- toJSON(out, auto_unbox = TRUE, digits = 3, na = "null")
dir.create("output/figures/interactive", showWarnings = FALSE, recursive = TRUE)
writeLines(js, "output/figures/interactive/interactive_data_rivers.json")
cat(sprintf("wrote interactive_data_rivers.json | sig=%d | all river cells=%d | size=%.1f MB\n",
    length(recs), length(ctx$lon),
    file.info("output/figures/interactive/interactive_data_rivers.json")$size/1e6))

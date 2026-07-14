# F1_phenomenon.R  -- AQUATIC PAPER, Figure 1: the phenomenon
# Establishes, intuitively and in physical units, that ENSO leaves a coherent, sign-flipping
# thermal imprint on inland waters -- the motivation the footprint (F2) then quantifies.
# Detrended composite = mean(detrended temp over El Nino years) - mean(over neutral years),
# so it is the SAME interannual signal the footprint pipeline uses (not the legacy raw masks).
# 2x2: {lake surface, river} x {El Nino (ENSO+), La Nina (ENSO-)}. Robinson.

suppressMessages({library(ncdf4); library(terra); library(tidyterra); library(ggplot2)
  library(sf); library(patchwork); library(rnaturalearth)})
## --- locate repo root (reproducible: Rscript, or source() in RStudio, from any dir) ---
local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  o <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  d <- if (length(a)) dirname(a[1]) else if (!is.null(o)) dirname(o) else getwd()
  d <- normalizePath(d, mustWork = FALSE)
  while (!all(dir.exists(file.path(d, c("input", "code")))) && dirname(d) != d) d <- dirname(d)
  assign("base", d, envir = .GlobalEnv)
})
lakes_dir <- base; rivers_dir <- base
setwd(base)
robin <- "+proj=robin +datum=WGS84 +no_defs"; crop_ll <- ext(-180,180,-58,85)
years <- 1901:2021; keep <- years>=1950 & years<=2021; yrs <- years[keep]; n <- length(yrs)

# linear-detrend operator (interannual band), identical to the footprint engine
X_t <- cbind(1, scale(yrs, scale=FALSE)); Hcomp <- diag(n) - X_t%*%solve(crossprod(X_t))%*%t(X_t)

# ENSO phase years from the z-scored DJF ONI (top / bottom quartile vs neutral middle)
si <- read.csv("input/clim_index/seasonal_indices.csv"); si <- si[match(yrs,si$year),]
e  <- si$enso_djf
pos <- e >= quantile(e,.75); neg <- e <= quantile(e,.25); neu <- !pos & !neg

# build a detrended composite raster (deg C anomaly) for one temperature field, masked to
# the cells where ENSO is detrended-significant (enso_sig from that system's engine)
composite <- function(nc_path, varname, sig_path){
  nc <- nc_open(nc_path); v <- ncvar_get(nc, varname); nc_close(nc); v <- v[,,keep]
  r  <- rast(lapply(seq_len(n), function(k){x<-rast(t(v[,,k]));ext(x)<-c(-180,180,-90,90);crs(x)<-"EPSG:4326";x}))
  tab<- as.data.frame(r, xy=TRUE, na.rm=TRUE); xy<-tab[,1:2]; cells<-cellFromXY(r[[1]],as.matrix(xy))
  Yd <- t(Hcomp %*% t(as.matrix(tab[,-(1:2)])))                 # detrended, deg C
  cp <- rowMeans(Yd[,pos]) - rowMeans(Yd[,neu])                 # El Nino composite
  cn <- rowMeans(Yd[,neg]) - rowMeans(Yd[,neu])                 # La Nina composite
  sig<- rast(sig_path)[cells][,1]; sigm <- !is.na(sig) & sig==1
  mk <- function(vals, only_sig){ o<-r[[1]]; values(o)<-NA
    idx<- if(only_sig) which(sigm) else seq_along(vals); o[cells[idx]]<-vals[idx]; crop(o,crop_ll)}
  list(pos_all=mk(cp,FALSE), pos_sig=mk(cp,TRUE), neg_all=mk(cn,FALSE), neg_sig=mk(cn,TRUE))
}

Lk <- composite(file.path(lakes_dir,"input/laketemp/gotm_20crv3-era5_obsclim_histsoc_default_surftemp_global_annual_1901_2021.nc"),
                "surftemp", file.path(lakes_dir,"output/footprint/surftemp/enso_sig.tiff"))
Rv <- composite(file.path(rivers_dir,"input/rivertemp/watergap2-2e_20crv3-era5_obsclim_histsoc_default_rivtemp_global_annual_1901_2021.nc"),
                "rivtemp",  file.path(rivers_dir,"output/footprint/rivtemp/enso_sig.tiff"))

LIM <- 0.8   # symmetric colour limit (deg C); anomalies capped for display
land <- ne_countries(scale="medium",returnclass="sf")|>st_make_valid()|>
  st_crop(xmin=-180,xmax=180,ymin=-58,ymax=85)|>st_transform(robin)
grat <- st_graticule(lat=seq(-30,60,30),lon=seq(-120,120,60))

panel <- function(all_r, sig_r, tag, title){
  cap <- function(x) clamp(x, -LIM, LIM, values=TRUE)
  ggplot()+
    geom_sf(data=land, fill="grey94", colour="grey86", linewidth=.1)+
    geom_spatraster(data=cap(all_r), alpha=.30)+
    geom_spatraster(data=cap(sig_r), alpha=1)+
    scale_fill_gradient2(name="Temperature anomaly (°C)", low="#2166AC", mid="grey96", high="#B2182B",
      midpoint=0, limits=c(-LIM,LIM), breaks=c(-LIM,0,LIM), labels=c("−0.8","0","+0.8"), na.value=NA)+
    geom_sf(data=grat, colour="grey88", linewidth=.15)+
    coord_sf(crs=robin, expand=FALSE)+
    labs(title=bquote(bold(.(tag))~"  "*.(title)))+
    theme_void(base_size=12)+
    theme(plot.title=element_text(size=12),legend.position="bottom",
      legend.key.width=unit(30,"pt"),legend.key.height=unit(9,"pt"),
      legend.title=element_text(size=9), plot.margin=margin(3,3,3,3))
}

fig <- (panel(Rv$pos_all,Rv$pos_sig,"a","Rivers · El Niño (ENSO+)") +
        panel(Rv$neg_all,Rv$neg_sig,"b","Rivers · La Niña (ENSO−)") +
        panel(Lk$pos_all,Lk$pos_sig,"c","Lakes · El Niño (ENSO+)") +
        panel(Lk$neg_all,Lk$neg_sig,"d","Lakes · La Niña (ENSO−)")) +
       plot_layout(ncol=2, guides="collect") & theme(legend.position="bottom")

dir.create("output/figures", showWarnings=FALSE, recursive=TRUE)
ggsave("output/figures/F1_phenomenon.png", fig, width=11, height=6.4, dpi=300, bg="white")
cat("saved F1_phenomenon.png\n")
cat(sprintf("El Nino years n=%d | La Nina n=%d | neutral n=%d\n", sum(pos),sum(neg),sum(neu)))

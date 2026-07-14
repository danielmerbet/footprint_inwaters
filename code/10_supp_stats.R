# 10_supp_stats.R  -- REVISION S-B(+)/S-C(+)/S-G: three referee-facing checks
#  (i)  index cross-correlation matrix + VIF on the DETRENDED, seasonal indices
#       (verify the near-1 VIF and expose any ENSO-IOD collinearity behind attribution)
#  (ii) frequency of NEGATIVE variation-partitioning components (adj-R^2 partitioning
#       can yield negative shared/unique fractions; 04 clamps medfrac -> we report the
#       raw incidence instead of hiding it)
#  (iii) interannual-band sensitivity: median F under FIRST-DIFFERENCING vs linear detrend
#       (a trend-free isolation of the interannual band, complementing the quadratic test)
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
X_t <- cbind(1, scale(yrs, scale=FALSE)); Hcomp <- diag(n) - X_t %*% solve(crossprod(X_t)) %*% t(X_t)
si <- read.csv("input/clim_index/seasonal_indices.csv"); si <- si[match(yrs, si$year), ]
outdir <- "output/figures/supplementary/robustness/"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
adjR2 <- function(y,X){p<-ncol(X);nn<-length(y);f<-.lm.fit(X,y);1-(sum(f$residuals^2)/sum((y-mean(y))^2))*(nn-1)/(nn-p)}

## (i) index collinearity -----------------------------------------------------
IDXraw <- cbind(ENSO=si$enso_djf, NAO=si$nao_djfm, IOD=si$iod_son)
IDXd   <- scale(Hcomp %*% IDXraw)
cormat <- cor(IDXd)
vif <- sapply(1:3, function(j) 1/(1-summary(lm(IDXd[,j]~IDXd[,-j]))$r.squared))
cat("=== (i) detrended seasonal index cross-correlation ===\n"); print(round(cormat,3))
cat("VIF:", paste(sprintf("%s=%.2f", colnames(IDXraw), vif), collapse="  "), "\n")
write.csv(round(cormat,3), paste0(outdir,"index_correlation.csv"))

## (ii) negative variation-partitioning components ---------------------------
cat("\n=== (ii) incidence of negative variation-partitioning components (sig cells) ===\n")
negtab <- list()
for (cp in c("surftemp","bottemp","rivtemp")) {
  d <- paste0("output/footprint/",cp,"/mechanism/")
  direct <- values(rast(paste0(d,"direct_adjR2.tiff")))[,1]
  medi   <- values(rast(paste0(d,"mediated_adjR2.tiff")))[,1]
  loc    <- values(rast(paste0(d,"local_adjR2.tiff")))[,1]
  ok <- !is.na(direct)
  r <- data.frame(compartment=cp, n=sum(ok),
                  pct_direct_neg = 100*mean(direct[ok] < 0),
                  pct_mediated_neg = 100*mean(medi[ok] < 0),
                  pct_local_neg = 100*mean(loc[ok] < 0))
  negtab[[cp]] <- r
  cat(sprintf("  %-9s n=%5d | direct<0: %.1f%% | mediated<0: %.1f%% | local<0: %.1f%%\n",
              cp, r$n, r$pct_direct_neg, r$pct_mediated_neg, r$pct_local_neg))
}
write.csv(do.call(rbind, negtab), paste0(outdir,"vp_negative_components.csv"), row.names=FALSE)

## (iii) first-difference interannual sensitivity ----------------------------
cat("\n=== (iii) median F: linear detrend vs first-difference (baseline-sig cells) ===\n")
tfiles <- list(
  surftemp=c("input/laketemp/gotm_20crv3-era5_obsclim_histsoc_default_surftemp_global_annual_1901_2021.nc","surftemp"),
  bottemp =c("input/laketemp/gotm_20crv3-era5_obsclim_histsoc_default_bottemp_global_annual_1901_2021.nc","bottemp"),
  rivtemp =c("input/rivertemp/watergap2-2e_20crv3-era5_obsclim_histsoc_default_rivtemp_global_annual_1901_2021.nc","rivtemp"))
IDXlin  <- IDXd
IDXdiff <- scale(apply(IDXraw, 2, diff))
sens <- list()
for (cp in names(tfiles)) {
  nc<-nc_open(tfiles[[cp]][1]); v<-ncvar_get(nc, tfiles[[cp]][2]); nc_close(nc); v<-v[,,keep]
  grid<-rast(lapply(seq_len(n),function(k){r<-rast(t(v[,,k]));ext(r)<-c(-180,180,-90,90);crs(r)<-"EPSG:4326";r}))
  tab<-as.data.frame(grid,xy=TRUE,na.rm=TRUE); cells<-cellFromXY(grid[[1]],as.matrix(tab[,1:2]))
  Yraw<-as.matrix(tab[,-(1:2)])
  q<-rast(paste0("output/footprint/",cp,"/qvalue_fdr.tiff"))[cells][,1]; sig<-which(q<0.05 & !is.na(q))
  Yd  <- t(Hcomp %*% t(Yraw)); Ydiff <- t(apply(Yraw,1,diff))
  mlin  <- median(sapply(sig, function(i) adjR2(Yd[i,],   cbind(1,IDXlin))),  na.rm=TRUE)
  mdiff <- median(sapply(sig, function(i) adjR2(Ydiff[i,],cbind(1,IDXdiff))), na.rm=TRUE)
  sens[[cp]] <- data.frame(compartment=cp, n_sig=length(sig), medF_linear=mlin, medF_firstdiff=mdiff)
  cat(sprintf("  %-9s n=%5d | linear detrend=%.3f | first-difference=%.3f\n", cp, length(sig), mlin, mdiff))
}
write.csv(do.call(rbind, sens), paste0(outdir,"firstdiff_sensitivity.csv"), row.names=FALSE)
cat("\nsaved -> ", outdir, "{index_correlation,vp_negative_components,firstdiff_sensitivity}.csv\n", sep="")

# 05_robustness.R  -- supplementary validation of the footprint analysis
#  S-A  block cross-validation (does F generalise, or overfit?)
#  S-B  collinearity / VIF for the index block and the local block
#  S-C  sensitivity of median footprint to lag scheme and detrend order
# Runs for lake surface and river (bottom analogous). Outputs supp figures + tables.

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
years <- 1901:2021; keep <- years >= 1950 & years <= 2021; yrs <- years[keep]; n <- length(yrs)
Hlin  <- function(){X<-cbind(1,poly(yrs,1));diag(n)-X%*%solve(crossprod(X))%*%t(X)}
Hquad <- function(){X<-cbind(1,poly(yrs,2));diag(n)-X%*%solve(crossprod(X))%*%t(X)}
detr  <- function(M,H) t(H%*%t(M))
si <- read.csv("input/clim_index/seasonal_indices.csv"); si <- si[match(yrs,si$year),]
mkidx <- function(H) scale(H%*%cbind(si$enso_djf,si$nao_djfm,si$iod_son))
adjR2 <- function(y,X){p<-ncol(X);nn<-length(y);f<-.lm.fit(X,y);1-(sum(f$residuals^2)/sum((y-mean(y))^2))*(nn-1)/(nn-p)}
vif   <- function(M){sapply(seq_len(ncol(M)),function(j){r2<-summary(lm(M[,j]~M[,-j]))$r.squared;1/(1-r2)})}

systems <- list(
  lakes  = list(nc="input/laketemp/gotm_20crv3-era5_obsclim_histsoc_default_surftemp_global_annual_1901_2021.nc", var="surftemp", fp="output/footprint/surftemp/", locv=c("tas","rsds","rlds","sfcwind"),      label="lake surface"),
  rivers = list(nc="input/rivertemp/watergap2-2e_20crv3-era5_obsclim_histsoc_default_rivtemp_global_annual_1901_2021.nc", var="rivtemp",  fp="output/footprint/rivtemp/",  locv=c("tas","rsds","rlds","sfcwind","pr"), label="river"))

outdir <- "output/figures/supplementary/robustness/"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
run_rob <- function(nm, cfg) {
  cat("\n########## robustness:", cfg$label, "##########\n")
  nc<-nc_open(cfg$nc); v<-ncvar_get(nc,cfg$var); nc_close(nc); v<-v[,,keep]
  grid<-rast(lapply(seq_len(n),function(k){r<-rast(t(v[,,k]));ext(r)<-c(-180,180,-90,90);crs(r)<-"EPSG:4326";r}))
  tab<-as.data.frame(grid,xy=TRUE,na.rm=TRUE); cells<-cellFromXY(grid[[1]],as.matrix(tab[,1:2])); Yraw<-as.matrix(tab[,-(1:2)])
  gv<-function(f) rast(paste0(cfg$fp,f))[cells][,1]
  q<-gv("qvalue_fdr.tiff"); lags<-cbind(gv("enso_lag.tiff"),gv("nao_lag.tiff"),gv("iod_lag.tiff")); sigr<-which(q<0.05 & !is.na(q))
  Hl<-Hlin(); Yd<-detr(Yraw,Hl); IDX<-mkidx(Hl); E<-IDX[,1];N<-IDX[,2];I<-IDX[,3]; fold<-cut(seq_len(n),6,labels=FALSE)
  cvF<-t(sapply(sigr,function(i){L<-lags[i,];mL<-max(L);rows<-(1+mL):n
    X<-cbind(1,E[rows-L[1]],N[rows-L[2]],I[rows-L[3]]); y<-Yd[i,rows]; fr<-fold[rows]; pred<-rep(NA,length(y))
    for(k in 1:6){te<-fr==k; if(sum(te)<2||sum(!te)<6)next; b<-.lm.fit(X[!te,,drop=FALSE],y[!te])$coefficients; pred[te]<-X[te,,drop=FALSE]%*%b}
    ok<-!is.na(pred); c(F=adjR2(y,X), cv=1-sum((y[ok]-pred[ok])^2)/sum((y[ok]-mean(y[ok]))^2))}))
  cat(sprintf("  S-A CV (n=%d): median F=%.3f | median CV R2=%.3f | %% CV>0 = %.1f\n",
      length(sigr),median(cvF[,"F"]),median(cvF[,"cv"]),100*mean(cvF[,"cv"]>0)))
  loadmet<-function(var){nc<-nc_open(paste0("input/clim_data/20crv3-era5_obsclim_",var,"_global_yearly_1901_2021.nc"))
    m<-ncvar_get(nc,var);nc_close(nc);m<-m[,,keep]
    st<-rast(lapply(seq_len(n),function(k){r<-rast(t(m[,,k]));ext(r)<-c(-180,180,-90,90);crs(r)<-"EPSG:4326";r}));detr(st[cells],Hl)}
  Lc<-lapply(cfg$locv,loadmet)
  vif_loc<-sapply(sigr,function(i){vif(sapply(Lc,function(x)x[i,]))}); rownames(vif_loc)<-cfg$locv
  cat("  S-B VIF indices:",paste(round(vif(IDX),2),collapse=","),"| local:",
      paste(sprintf("%s=%.1f",cfg$locv,apply(vif_loc,1,median)),collapse=" "),"\n")
  medF<-function(H){Yd2<-detr(Yraw,H);IDX2<-mkidx(H);sapply(sigr,function(i){L<-lags[i,];mL<-max(L);rows<-(1+mL):n
    adjR2(Yd2[i,rows],cbind(1,IDX2[rows-L[1],1],IDX2[rows-L[2],2],IDX2[rows-L[3],3]))})}
  medF0<-function(H){Yd2<-detr(Yraw,H);IDX2<-mkidx(H);sapply(sigr,function(i) adjR2(Yd2[i,],cbind(1,IDX2[,1],IDX2[,2],IDX2[,3])))}
  sens<-data.frame(scheme=c("baseline","zero lag","quadratic detrend"),median_F=c(median(medF(Hl)),median(medF0(Hl)),median(medF(Hquad()))))
  cat("  S-C sensitivity:",paste(sprintf("%.3f",sens$median_F),collapse=" / "),"\n")
  write.csv(sens,paste0(outdir,"sensitivity_",nm,".csv"),row.names=FALSE)
  p1<-ggplot(data.frame(F=cvF[,"F"],cv=cvF[,"cv"]),aes(F,cv))+geom_hline(yintercept=0,colour="grey80")+
    geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey60")+geom_point(alpha=.15,size=.7,colour="#21918c")+
    coord_equal(xlim=c(0,.7),ylim=c(-.2,.7))+labs(title=paste0("a  Cross-validation (",cfg$label,")"),
    x="in-sample F (adj. R²)",y="block-CV R²")+theme_minimal(base_size=11)
  p2<-ggplot(data.frame(var=factor(cfg$locv,levels=cfg$locv),vif=apply(vif_loc,1,median)),aes(var,vif))+
    geom_col(fill="#3b528b",width=.6)+geom_hline(yintercept=5,linetype="dashed",colour="#C1272D")+
    labs(title="b  Local-block VIF",x=NULL,y="median VIF")+theme_minimal(base_size=11)
  ggsave(paste0(outdir,"robustness_",nm,".png"),p1+p2+plot_layout(widths=c(1.2,1)),width=10,height=4.4,dpi=200,bg="white")
}
for(nm in names(systems)) run_rob(nm, systems[[nm]])
cat("\nsaved ->",outdir,"\n")

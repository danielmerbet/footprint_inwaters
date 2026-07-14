# FS_nao_iod.R  -- SUPPLEMENTARY: NAO and IOD analogues of Figure 1 (phenomenon) and of
# Figure 4b (memory axis) + 4d (contemporaneous response). ENSO is the clearest mode but
# NAO and IOD organise important regional footprints; this script documents them with the
# SAME diagnostics used for ENSO in the main text.
#   Outputs (output/figures/supplementary/):
#     FS_nao_phenomenon.png, FS_iod_phenomenon.png   (2x2 composite maps, sign-flip)
#     FS_nao_synthesis.png,  FS_iod_synthesis.png    (memory axis + lag-0 response)
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
setwd(base)
outdir <- "output/figures/supplementary/"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
robin <- "+proj=robin +datum=WGS84 +no_defs"; crop_ll <- ext(-180,180,-58,85)
years <- 1901:2021; keep <- years>=1950 & years<=2021; yrs <- years[keep]; n <- length(yrs)
X_t <- cbind(1, scale(yrs, scale=FALSE)); Hcomp <- diag(n) - X_t%*%solve(crossprod(X_t))%*%t(X_t)
si <- read.csv("input/clim_index/seasonal_indices.csv"); si <- si[match(yrs,si$year),]

surf_nc <- "input/laketemp/gotm_20crv3-era5_obsclim_histsoc_default_surftemp_global_annual_1901_2021.nc"
bot_nc  <- "input/laketemp/gotm_20crv3-era5_obsclim_histsoc_default_bottemp_global_annual_1901_2021.nc"
riv_nc  <- "input/rivertemp/watergap2-2e_20crv3-era5_obsclim_histsoc_default_rivtemp_global_annual_1901_2021.nc"

load_detr <- function(nc_path, varname) {                       # -> list(Yd, xy, cells, grid)
  nc <- nc_open(nc_path); v <- ncvar_get(nc, varname); nc_close(nc); v <- v[,,keep]
  r  <- rast(lapply(seq_len(n), function(k){x<-rast(t(v[,,k]));ext(x)<-c(-180,180,-90,90);crs(x)<-"EPSG:4326";x}))
  tab<- as.data.frame(r, xy=TRUE, na.rm=TRUE)
  list(Yd=t(Hcomp%*%t(as.matrix(tab[,-(1:2)]))), xy=tab[,1:2],
       cells=cellFromXY(r[[1]],as.matrix(tab[,1:2])), grid=r)
}

##====================== PART 1: phenomenon composites =========================
land <- ne_countries(scale="medium",returnclass="sf")|>st_make_valid()|>
  st_crop(xmin=-180,xmax=180,ymin=-58,ymax=85)|>st_transform(robin)
grat <- st_graticule(lat=seq(-30,60,30),lon=seq(-120,120,60))
LIM  <- 0.8
panel_map <- function(all_r, sig_r, tag, title){
  cap <- function(x) clamp(x, -LIM, LIM, values=TRUE)
  ggplot()+ geom_sf(data=land, fill="grey94", colour="grey86", linewidth=.1)+
    geom_spatraster(data=cap(all_r), alpha=.30)+ geom_spatraster(data=cap(sig_r), alpha=1)+
    scale_fill_gradient2(name="Temperature anomaly (°C)", low="#2166AC", mid="grey96", high="#B2182B",
      midpoint=0, limits=c(-LIM,LIM), breaks=c(-LIM,0,LIM), labels=c("−0.8","0","+0.8"), na.value=NA)+
    geom_sf(data=grat, colour="grey88", linewidth=.15)+ coord_sf(crs=robin, expand=FALSE)+
    labs(title=bquote(bold(.(tag))~"  "*.(title)))+ theme_void(base_size=12)+
    theme(plot.title=element_text(size=12),legend.position="bottom",
      legend.key.width=unit(30,"pt"),legend.key.height=unit(9,"pt"),
      legend.title=element_text(size=9), plot.margin=margin(3,3,3,3))
}
composite <- function(dd, idx, sig_path){
  pos <- idx>=quantile(idx,.75); neg <- idx<=quantile(idx,.25); neu <- !pos & !neg
  cp <- rowMeans(dd$Yd[,pos]) - rowMeans(dd$Yd[,neu])
  cn <- rowMeans(dd$Yd[,neg]) - rowMeans(dd$Yd[,neu])
  sig<- rast(sig_path)[dd$cells][,1]; sigm <- !is.na(sig) & sig==1
  mk <- function(vals, only_sig){ o<-dd$grid[[1]]; values(o)<-NA
    idx2<- if(only_sig) which(sigm) else seq_along(vals); o[dd$cells[idx2]]<-vals[idx2]; crop(o,crop_ll)}
  list(pos_all=mk(cp,FALSE),pos_sig=mk(cp,TRUE),neg_all=mk(cn,FALSE),neg_sig=mk(cn,TRUE),
       npos=sum(pos),nneg=sum(neg),nneu=sum(neu))
}
make_phenomenon <- function(idxcol, sfx, name, plus, minus, fout){
  dS <- load_detr(surf_nc,"surftemp"); dR <- load_detr(riv_nc,"rivtemp")
  Lk <- composite(dS, si[[idxcol]], sprintf("output/footprint/surftemp/%s_sig.tiff", sfx))
  Rv <- composite(dR, si[[idxcol]], sprintf("output/footprint/rivtemp/%s_sig.tiff",  sfx))
  fig <- (panel_map(Rv$pos_all,Rv$pos_sig,"a",sprintf("Rivers · %s",plus)) +
          panel_map(Rv$neg_all,Rv$neg_sig,"b",sprintf("Rivers · %s",minus)) +
          panel_map(Lk$pos_all,Lk$pos_sig,"c",sprintf("Lakes · %s",plus)) +
          panel_map(Lk$neg_all,Lk$neg_sig,"d",sprintf("Lakes · %s",minus))) +
         plot_layout(ncol=2, guides="collect") & theme(legend.position="bottom")
  ggsave(fout, fig, width=11, height=6.4, dpi=300, bg="white")
  cat(sprintf("saved %s | %s+ n=%d  %s- n=%d  neutral n=%d\n", basename(fout), name,Lk$npos,name,Lk$nneg,Lk$nneu))
}
make_phenomenon("nao_djfm","nao","NAO","NAO+ (positive)","NAO− (negative)", paste0(outdir,"FS_nao_phenomenon.png"))
make_phenomenon("iod_son","iod","IOD","IOD+ (positive)","IOD− (negative)", paste0(outdir,"FS_iod_phenomenon.png"))

##====================== PART 2: memory axis + lag-0 response ==================
SYS <- list(
  list(key="River",        dir="output/footprint/rivtemp/",  nc=riv_nc,  var="rivtemp",  col="#E08A00"),
  list(key="Lake surface", dir="output/footprint/surftemp/", nc=surf_nc, var="surftemp", col="#4A90B8"),
  list(key="Lake bottom",  dir="output/footprint/bottemp/",  nc=bot_nc,  var="bottemp",  col="#16324F"))
lev <- sapply(SYS,`[[`,"key"); pal <- setNames(sapply(SYS,`[[`,"col"),lev); mem_order <- c("River","Lake surface","Lake bottom")
g <- function(s,f) values(rast(paste0(s$dir,f)))[,1]

# target: which dominant-index column (ENSO=1,IOD=2,NAO=3) and its detrended standardized series
synth_for <- function(target_col, iseries, iname, ilabel, fout){
  ISER <- as.numeric(scale(Hcomp %*% iseries))
  per <- lapply(SYS, function(s){
    q<-g(s,"qvalue_fdr.tiff"); sig<-which(q<0.05 & !is.na(q))
    E<-g(s,"enso_adjR2.tiff"); I<-g(s,"iod_adjR2.tiff"); N<-g(s,"nao_adjR2.tiff")
    dom<-max.col(cbind(E,I,N), ties.method="first")           # 1=ENSO,2=IOD,3=NAO
    mf<-g(s,"mechanism/mediated_frac.tiff")
    domsig <- sig[dom[sig]==target_col]
    direct <- pmax(0, 1 - mf[domsig])
    dd<-load_detr(s$nc,s$var); ord<-match(domsig, dd$cells); ord<-ord[!is.na(ord)]
    Yd<-dd$Yd[ord,,drop=FALSE]
    sgn<-sign(apply(Yd,1,function(y) suppressWarnings(cor(y,ISER)))); sgn[is.na(sgn)]<-1
    comp<-as.numeric(scale(colMeans(Yd*sgn)))                  # sign-aligned congruent composite
    lags<- -3:3; cc<-sapply(lags,function(L){a<-comp;b<-ISER
      if(L>0){a<-a[(1+L):n];b<-b[1:(n-L)]} else if(L<0){a<-a[1:(n+L)];b<-b[(1-L):n]}; cor(a,b)})
    list(key=s$key, direct=direct, comp=comp, r0=cor(comp,ISER), peaklag=lags[which.max(abs(cc))], ndom=length(domsig))
  })
  names(per)<-lev
  cat(sprintf("\n[%s] dominant-cell response (sign-aligned):\n", iname))
  for(p in per) cat(sprintf("   %-12s n=%4d | r(lag0)=%.2f | CCF peak lag=%d | median direct=%.3f\n",
                            p$key,p$ndom,p$r0,p$peaklag, median(p$direct)))
  th<-theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),
    plot.title=element_text(size=12,face="bold"),legend.position="none")
  # (a) memory axis among this index's dominant cells
  dfl<-do.call(rbind,lapply(per,function(x) data.frame(system=x$key,direct=x$direct)))
  dfl$system<-factor(dfl$system,levels=mem_order)
  pb<-ggplot(dfl,aes(system,direct,fill=system))+
    geom_violin(colour=NA,alpha=.55,scale="width")+
    geom_boxplot(width=.13,outlier.shape=NA,colour="grey25",fill="white",linewidth=.4)+
    scale_fill_manual(values=pal)+coord_cartesian(ylim=c(0,0.6))+
    labs(title=sprintf("a  Memory axis · %s-dominated cells",iname),x=NULL,y="Direct (memory) fraction of F")+
    th+theme(axis.text.x=element_text(angle=15,hjust=1))
  # (b) contemporaneous response
  tsl<-do.call(rbind,lapply(per,function(x) data.frame(year=yrs,system=x$key,z=x$comp)))
  tsl$system<-factor(tsl$system,levels=lev); ens<-data.frame(year=yrs,z=ISER)
  pd<-ggplot()+ geom_col(data=ens,aes(year,z),fill="grey85",width=.85)+
    geom_line(data=tsl,aes(year,z,colour=system),linewidth=.7)+ scale_colour_manual(values=pal,name=NULL)+
    labs(title=sprintf("b  Contemporaneous response (lag 0) · %s",iname),x=NULL,y="Standardized anomaly (index-congruent)")+
    annotate("text",x=min(yrs)+3,y=2.6,label=sprintf("grey bars = %s",ilabel),size=2.9,colour="grey45",hjust=0)+
    th+theme(legend.position="top",legend.text=element_text(size=8.5),legend.key.size=unit(12,"pt"))
  ggsave(fout, pb|pd, width=10, height=4.3, dpi=300, bg="white")
  cat(sprintf("saved %s\n", basename(fout)))
}
synth_for(3, si$nao_djfm, "NAO", "NAO (DJFM)", paste0(outdir,"FS_nao_synthesis.png"))
synth_for(2, si$iod_son,  "IOD", "IOD (SON)",  paste0(outdir,"FS_iod_synthesis.png"))
cat("\ndone -> ", outdir, "FS_{nao,iod}_{phenomenon,synthesis}.png\n", sep="")

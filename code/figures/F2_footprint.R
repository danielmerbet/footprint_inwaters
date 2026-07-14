# F2_footprint.R  -- AQUATIC PAPER, Figure 2: the quantified footprint
# 3 systems (lake surface, lake bottom, river) x 2 columns:
#   left  = oscillation footprint F (adj R^2), FDR-significant sentinels, viridis
#   right = dominant oscillation (ENSO/IOD/NAO) among sentinels, Okabe-Ito
# This QUANTIFIES the sign-pattern F1 established. Merges the old fig3 + fig3b.

suppressMessages({library(terra); library(tidyterra); library(ggplot2); library(sf)
  library(viridis); library(patchwork); library(rnaturalearth)})
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
robin <- "+proj=robin +datum=WGS84 +no_defs"; crop_ll <- ext(-180,180,-58,85); FMAX <- 0.6
cols <- c(ENSO="#0072B2", IOD="#D55E00", NAO="#009E73")
land <- ne_countries(scale="medium",returnclass="sf")|>st_make_valid()|>
  st_crop(xmin=-180,xmax=180,ymin=-58,ymax=85)|>st_transform(robin)
grat <- st_graticule(lat=seq(-30,60,30),lon=seq(-120,120,60))
base_gg <- function() list(geom_sf(data=land,fill="grey94",colour="grey86",linewidth=.1),
  geom_sf(data=grat,colour="grey88",linewidth=.15), coord_sf(crs=robin,expand=FALSE),
  theme_void(base_size=12),
  theme(plot.title=element_text(size=12),plot.subtitle=element_text(size=9.5,colour="grey30"),
    legend.position="bottom",legend.title=element_text(size=9),
    legend.key.width=unit(34,"pt"),legend.key.height=unit(9,"pt"),plot.margin=margin(3,4,3,4)))

Fpanel <- function(d,tag,title){
  F<-crop(rast(paste0(d,"footprint_adjR2.tiff")),crop_ll); q<-crop(rast(paste0(d,"qvalue_fdr.tiff")),crop_ll)
  Fsig<-F; Fsig[q>=0.05|is.na(q)]<-NA; Fcap<-clamp(Fsig,upper=FMAX,values=TRUE); Fbg<-ifel(!is.na(F),1,NA)
  medF<-round(median(values(Fsig),na.rm=TRUE),2); pct<-round(100*sum(!is.na(values(Fsig)))/sum(!is.na(values(F))),1)
  ggplot()+base_gg()[1]+
    geom_spatraster(data=Fbg)+scale_fill_gradient(low="grey87",high="grey87",na.value=NA,guide="none")+
    ggnewscale::new_scale_fill()+ geom_spatraster(data=Fcap)+
    scale_fill_viridis_c(name="Footprint  F  (adj. R²)",limits=c(0,FMAX),breaks=c(0,0.3,0.6),na.value=NA)+
    base_gg()[-1]+ labs(title=bquote(bold(.(tag))~"  "*.(title)),
      subtitle=bquote("median F = "*.(medF)*"  ·  "*.(pct)*"% teleconnected"))
}
Dpanel <- function(d,tag,title){
  q<-crop(rast(paste0(d,"qvalue_fdr.tiff")),crop_ll)
  st<-c(crop(rast(paste0(d,"enso_adjR2.tiff")),crop_ll),crop(rast(paste0(d,"iod_adjR2.tiff")),crop_ll),
        crop(rast(paste0(d,"nao_adjR2.tiff")),crop_ll))
  dom<-which.max(st); dom[q>=0.05|is.na(q)]<-NA; levels(dom)<-data.frame(id=1:3,index=c("ENSO","IOD","NAO"))
  Fbg<-ifel(!is.na(q),1,NA); tb<-round(100*table(factor(values(dom),1:3))/sum(!is.na(values(dom))))
  ggplot()+base_gg()[1]+
    geom_spatraster(data=Fbg)+scale_fill_gradient(low="grey87",high="grey87",na.value=NA,guide="none")+
    ggnewscale::new_scale_fill()+ geom_spatraster(data=dom)+
    scale_fill_manual(values=cols,na.value=NA,na.translate=FALSE,name="Dominant oscillation",drop=FALSE)+
    base_gg()[-1]+ labs(title=bquote(bold(.(tag))~"  "*.(title)),
      subtitle=sprintf("ENSO %d%% · IOD %d%% · NAO %d%%",tb[1],tb[2],tb[3]))
}
FP_OUT_ROOT <- Sys.getenv("FP_OUT_ROOT","output/footprint/")   # override for annual-index supp. variant
dl<-paste0(lakes_dir,"/",FP_OUT_ROOT,"surftemp/")
db<-paste0(lakes_dir,"/",FP_OUT_ROOT,"bottemp/")
dr<-paste0(rivers_dir,"/",FP_OUT_ROOT,"rivtemp/")

fig <- (Fpanel(dr,"a","River")        | Dpanel(dr,"b","River")) /
       (Fpanel(dl,"c","Lake surface") | Dpanel(dl,"d","Lake surface")) /
       (Fpanel(db,"e","Lake bottom")  | Dpanel(db,"f","Lake bottom"))  +
       plot_layout(guides="collect") & theme(legend.position="bottom")

FIG_OUT <- Sys.getenv("FIG_OUT","output/figures/F2_footprint.png")
dir.create(dirname(FIG_OUT),showWarnings=FALSE,recursive=TRUE)
ggsave(FIG_OUT,fig,width=11,height=9.6,dpi=300,bg="white")
cat("saved",FIG_OUT,"\n")

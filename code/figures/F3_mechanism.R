# F3_mechanism.R  -- AQUATIC PAPER, Figure 3: how the signal arrives (REPLACES old CatBoost fig2)
# Per sentinel cell, the dominant pathway by which the oscillation reaches the water:
#   air temperature / solar / longwave / wind mixing / precipitation (rivers) / direct (memory).
# 3 stacked panels: lake surface, lake bottom, river. Unified 6-category scheme.
# Two-block variation partitioning (mediation_partitioning.R); "direct" where mediated frac < 0.7.

suppressMessages({library(terra); library(tidyterra); library(ggplot2); library(sf); library(patchwork); library(rnaturalearth)})
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
robin<-"+proj=robin +datum=WGS84 +no_defs"; crop_ll<-ext(-180,180,-58,85); DIRECT_T<-0.7
labs6<-c("Air temperature","Solar radiation","Longwave radiation","Wind mixing","Precipitation","Direct (memory)")
cols6<-c("Air temperature"="#C1272D","Solar radiation"="#E8A33D","Longwave radiation"="#6A4C9C",
         "Wind mixing"="#2E86AB","Precipitation"="#1B9E77","Direct (memory)"="#1A1A1A")
land<-ne_countries(scale="medium",returnclass="sf")|>st_make_valid()|>
  st_crop(xmin=-180,xmax=180,ymin=-58,ymax=85)|>st_transform(robin)
grat<-st_graticule(lat=seq(-30,60,30),lon=seq(-120,120,60))

# for lakes dominant_mediator in {1..4} (tas,rsds,rlds,sfcwind); rivers {1..5} (+pr).
# Unify: keep those ids, and set DIRECT to id 6 wherever mediated fraction < DIRECT_T.
panel<-function(vp_dir, tag, title){
  m<-paste0(vp_dir,"mechanism/")
  dm<-crop(rast(paste0(m,"dominant_mediator.tiff")),crop_ll)
  mf<-crop(rast(paste0(m,"mediated_frac.tiff")),crop_ll)
  cat_r<-dm; cat_r[mf<DIRECT_T]<-6
  levels(cat_r)<-data.frame(id=1:6,mech=labs6)
  tb<-round(100*table(factor(values(cat_r),1:6))/sum(!is.na(values(cat_r))))
  sub<-sprintf("air-temp %d%% · longwave %d%% · direct %d%%",tb[1],tb[3],tb[6])
  q<-crop(rast(paste0(vp_dir,"qvalue_fdr.tiff")),crop_ll); Fbg<-ifel(!is.na(q),1,NA)
  ggplot()+
    geom_sf(data=land,fill="grey93",colour="grey85",linewidth=.1)+
    geom_spatraster(data=Fbg)+scale_fill_gradient(low="grey86",high="grey86",na.value=NA,guide="none")+
    ggnewscale::new_scale_fill()+ geom_spatraster(data=cat_r)+
    scale_fill_manual(values=cols6,na.value=NA,na.translate=FALSE,limits=labs6,drop=FALSE,guide="none")+
    geom_sf(data=grat,colour="grey88",linewidth=.15)+ coord_sf(crs=robin,expand=FALSE)+
    labs(title=bquote(bold(.(tag))~"  "*.(title)), subtitle=sub)+
    theme_void(base_size=12)+
    theme(plot.title=element_text(size=12),plot.subtitle=element_text(size=9.5,colour="grey30"),
      legend.position="none",plot.margin=margin(3,4,3,4))
}
# standalone legend row (drawn explicitly to avoid patchwork guide-collection duplication)
legdf<-data.frame(lab=factor(labs6,levels=labs6), x=c(0,1.15,2.5,0,1.15,2.5), y=c(1,1,1,0,0,0))
pleg<-ggplot(legdf,aes(x,y))+
  geom_point(aes(colour=lab),size=6,shape=15)+
  geom_text(aes(label=lab),hjust=0,nudge_x=.09,size=3.5)+
  scale_colour_manual(values=cols6,guide="none")+
  annotate("text",x=0,y=1.75,label="How the signal arrives",fontface="bold",size=3.6,hjust=0)+
  scale_x_continuous(limits=c(-.1,3.7))+scale_y_continuous(limits=c(-.6,2.1))+
  theme_void()
fig<-( panel(paste0(rivers_dir,"/output/footprint/rivtemp/"), "a","River") /
       panel(paste0(lakes_dir,"/output/footprint/surftemp/"),"b","Lake surface") /
       panel(paste0(lakes_dir,"/output/footprint/bottemp/"), "c","Lake bottom")  /
       pleg ) + plot_layout(heights=c(1,1,1,0.24))

dir.create("output/figures",showWarnings=FALSE,recursive=TRUE)
ggsave("output/figures/F3_mechanism.png",fig,width=8.6,height=8.8,dpi=300,bg="white")
cat("saved F3_mechanism.png\n")

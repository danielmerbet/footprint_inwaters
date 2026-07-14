# F4_synthesis.R  -- AQUATIC PAPER, Figure 4: exposed surface waters vs the buffered deep lake
# The maps (F2/F3) imply a spectrum; this figure states it, on the axes the DATA supports.
# NOTE: the detrended engine assigns lag 0 to ~all sentinel cells (ENSO/NAO contemporaneous
# everywhere; IOD lag-1 in ~2%), so the contrast is NOT speed/lag -- it is MEMORY (direct
# fraction) and amplitude. Rivers & lake surfaces are exposed, fully-mediated integrators;
# only the deep lake develops a direct/memory pathway.
#   (a) strength vs extent
#   (b) memory axis        - direct-fraction (1 - mediated) distribution per compartment
#   (c) mediator identity  - which local pathway carries the signal (incl. river precipitation)
#   (d) contemporaneous response - tropical-ENSO composite series: all track ENSO at lag 0,
#                                  the deep lake damped (buffered), not delayed
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
lakes_dir <- base; rivers_dir <- base
setwd(base)
years<-1901:2021; keep<-years>=1950&years<=2021; yrs<-years[keep]; n<-length(yrs)
X_t<-cbind(1,scale(yrs,scale=FALSE)); Hcomp<-diag(n)-X_t%*%solve(crossprod(X_t))%*%t(X_t)
si<-read.csv("input/clim_index/seasonal_indices.csv"); si<-si[match(yrs,si$year),]
ENSO<-as.numeric(scale(Hcomp%*%si$enso_djf))

SYS<-list(
  list(key="River",dir=paste0(rivers_dir,"/output/footprint/rivtemp/"),col="#E08A00",
       nc=paste0(rivers_dir,"/input/rivertemp/watergap2-2e_20crv3-era5_obsclim_histsoc_default_rivtemp_global_annual_1901_2021.nc"),var="rivtemp"),
  list(key="Lake surface",dir=paste0(lakes_dir,"/output/footprint/surftemp/"),col="#4A90B8",
       nc=paste0(lakes_dir,"/input/laketemp/gotm_20crv3-era5_obsclim_histsoc_default_surftemp_global_annual_1901_2021.nc"),var="surftemp"),
  list(key="Lake bottom",dir=paste0(lakes_dir,"/output/footprint/bottemp/"),col="#16324F",
       nc=paste0(lakes_dir,"/input/laketemp/gotm_20crv3-era5_obsclim_histsoc_default_bottemp_global_annual_1901_2021.nc"),var="bottemp"))
lev<-sapply(SYS,`[[`,"key"); pal<-setNames(sapply(SYS,`[[`,"col"),lev)
medlab6<-c("Air temperature","Solar radiation","Longwave radiation","Wind mixing","Precipitation","Direct (memory)")
medcol6<-c("Air temperature"="#C1272D","Solar radiation"="#E8A33D","Longwave radiation"="#6A4C9C",
           "Wind mixing"="#2E86AB","Precipitation"="#1B9E77","Direct (memory)"="#1A1A1A")

grab<-function(s){
  r0<-rast(paste0(s$dir,"footprint_adjR2.tiff")); g<-function(f) values(rast(paste0(s$dir,f)))[,1]
  q<-g("qvalue_fdr.tiff"); F<-g("footprint_adjR2.tiff"); sig<-which(q<0.05&!is.na(q))
  E<-g("enso_adjR2.tiff"); I<-g("iod_adjR2.tiff"); N<-g("nao_adjR2.tiff")
  dom<-max.col(cbind(E,I,N),ties.method="first")
  mf<-values(rast(paste0(s$dir,"mechanism/mediated_frac.tiff")))[,1]
  dm<-values(rast(paste0(s$dir,"mechanism/dominant_mediator.tiff")))[,1]
  cat6<-dm; cat6[!is.na(mf)&mf<0.7]<-6                     # unify: direct -> 6
  # tropical ENSO-dominant teleconnected composite series (detrended, standardized)
  nc<-nc_open(s$nc); v<-ncvar_get(nc,s$var); nc_close(nc); v<-v[,,keep]
  rr<-rast(lapply(seq_len(n),function(k){x<-rast(t(v[,,k]));ext(x)<-c(-180,180,-90,90);crs(x)<-"EPSG:4326";x}))
  tab<-as.data.frame(rr,xy=TRUE,na.rm=TRUE); cells<-cellFromXY(rr[[1]],as.matrix(tab[,1:2]))
  ord<-match(sig,cells); trop<-which(abs(tab$y[ord])<23.5 & dom[sig]==1 & !is.na(ord))
  Yd<-t(Hcomp%*%t(as.matrix(tab[ord[trop],-(1:2)])))
  comp<-as.numeric(scale(colMeans(Yd)))
  list(key=s$key,pct=100*length(sig)/sum(!is.na(F)),medF=median(F[sig]),q25=quantile(F[sig],.25),q75=quantile(F[sig],.75),
       direct=1-mf[sig], cat6=cat6[sig], comp=comp, rENSO=cor(comp,ENSO), ntrop=length(trop))
}
D<-lapply(SYS,grab); names(D)<-lev
sm<-data.frame(system=factor(lev,levels=lev),pct=sapply(D,`[[`,"pct"),medF=sapply(D,`[[`,"medF"),
  q25=sapply(D,`[[`,"q25"),q75=sapply(D,`[[`,"q75"),rENSO=sapply(D,`[[`,"rENSO"),ntrop=sapply(D,`[[`,"ntrop"))
cat("summary:\n"); print(sm,row.names=FALSE)

th<-theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),
  plot.title=element_text(size=12,face="bold"),legend.position="none")

## (a) strength vs extent
pa<-ggplot(sm,aes(pct,medF,colour=system))+
  geom_errorbar(aes(ymin=q25,ymax=q75),width=.28,linewidth=.6)+ geom_point(size=4.2)+
  ggrepel::geom_text_repel(aes(label=system),size=3.4,fontface="bold",min.segment.length=99,point.padding=9)+
  scale_colour_manual(values=pal)+ scale_x_continuous(limits=c(0,10))+
  labs(title="a  Strength vs extent",x="Teleconnected cells (% of tested)",y="Footprint F  (median, IQR)")+th

## (b) memory axis: direct-fraction distribution (river, lake surface, lake bottom)
mem_order<-c("River","Lake surface","Lake bottom")
dfl<-do.call(rbind,lapply(D,function(x) data.frame(system=x$key,direct=pmax(0,x$direct))))
dfl$system<-factor(dfl$system,levels=mem_order)
pb<-ggplot(dfl,aes(system,direct,fill=system))+
  geom_violin(colour=NA,alpha=.55,scale="width")+
  geom_boxplot(width=.13,outlier.shape=NA,colour="grey25",fill="white",linewidth=.4)+
  scale_fill_manual(values=pal)+ coord_cartesian(ylim=c(0,0.6))+
  labs(title="b  Memory axis",x=NULL,y="Direct (memory) fraction of F")+
  annotate("text",x=3,y=.55,label="deep water\nstores memory",size=2.9,colour="grey35",lineheight=.9)+
  th+theme(axis.text.x=element_text(angle=15,hjust=1))

## (c) mediator identity composition
cml<-do.call(rbind,lapply(D,function(x){t<-table(factor(x$cat6,levels=1:6))
  data.frame(system=x$key,mech=factor(medlab6,levels=medlab6),frac=100*as.numeric(t)/sum(t))}))
cml$system<-factor(cml$system,levels=mem_order)
pc<-ggplot(cml,aes(system,frac,fill=mech))+geom_col(width=.72)+
  scale_fill_manual(values=medcol6,name=NULL)+
  labs(title="c  Mediator identity",x=NULL,y="Teleconnected cells (%)")+
  th+theme(legend.position="right",legend.text=element_text(size=8),legend.key.size=unit(9,"pt"),
    axis.text.x=element_text(angle=15,hjust=1))

## (d) contemporaneous tropical-ENSO composite response (lag 0; deep water damped)
tsl<-do.call(rbind,lapply(D,function(x) data.frame(year=yrs,system=x$key,z=x$comp)))
tsl$system<-factor(tsl$system,levels=lev)
ens<-data.frame(year=yrs,z=ENSO)
pd<-ggplot()+
  geom_col(data=ens,aes(year,z),fill="grey85",width=.85)+
  geom_line(data=tsl,aes(year,z,colour=system),linewidth=.7)+
  scale_colour_manual(values=pal,name=NULL)+
  labs(title="d  Contemporaneous response (lag 0)",
       x=NULL,y="Standardized anomaly")+
  annotate("text",x=1954,y=2.6,label="grey bars = ENSO (DJF)",size=2.9,colour="grey45",hjust=0)+
  th+theme(legend.position="top",legend.text=element_text(size=8.5),legend.key.size=unit(12,"pt"))

fig<-(pa|pb)/((pc|pd)+plot_layout(widths=c(1,1.5)))
dir.create("output/figures",showWarnings=FALSE,recursive=TRUE)
ggsave("output/figures/F4_synthesis.png",fig,width=10.4,height=8.6,dpi=300,bg="white")
cat("saved F4_synthesis.png\n")
cat(sprintf("tropical-ENSO composite corr with ENSO: %s\n",
  paste(sprintf("%s r=%.2f (n=%d)",lev,sm$rENSO,sm$ntrop),collapse=" | ")))

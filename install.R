# install.R — R dependencies for "The footprint of climate oscillations in global inland waters"
# Run once:  Rscript install.R   (or  source("install.R")  in an R session)
#
# System libraries required first (install via your OS package manager):
#   GDAL, GEOS, PROJ            -> for {terra}, {sf}
#   NetCDF C library (libnetcdf) -> for {ncdf4}
# Optional external tools:
#   cdo (Climate Data Operators) -> NetCDF compression during data preparation
#   pdflatex (TeX Live/MiKTeX)    -> to build the manuscript

pkgs <- c(
  "terra",         # raster I/O and grid algebra (GeoTIFF, NetCDF)
  "ncdf4",         # NetCDF reading
  "sf",            # vector geometry (coastlines, graticules)
  "ggplot2",       # figures
  "patchwork",     # figure composition
  "jsonlite",      # interactive-explorer JSON export
  "tidyterra",     # ggplot geoms for terra rasters
  "rnaturalearth", # basemap land polygons
  "viridis",       # perceptually-uniform colour scales
  "ggrepel",       # non-overlapping labels (F4)
  "ggnewscale"     # multiple fill scales in one map (F2/F3)
)

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All required packages already installed.")
}

# --- Tested environment (developed & verified on) -------------------------------
# R 4.1.2
#   terra 1.8.60 · ncdf4 1.23 · sf 1.0.16 · ggplot2 3.5.2 · patchwork 1.3.0
#   jsonlite 1.9.1 · tidyterra 0.7.2 · rnaturalearth 1.1.0
# Record your own environment for provenance with:
#   writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

# build_annual_rivtemp.R
# Aggregate WaterGAP2-2e monthly river temperature (triver) to annual means,
# writing an annual NetCDF that mirrors the lake annual file so the footprint
# engine runs unchanged. Variable renamed triver -> rivtemp.
# Memory-safe: reads the monthly file one year (12-month slab) at a time,
# never holding the full [720,360,1452] array in RAM.
suppressMessages(library(ncdf4))
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
# Source monthly file is kept out of the repo (too large even compressed, ~236MB);
# it lives in largefiles/. The annual output below is committed instead.
fin <- "largefiles/watergap2-2e_20crv3-era5_obsclim_histsoc_default_triver_global_monthly_1901_2021.nc"
nc  <- nc_open(fin)
lon <- ncvar_get(nc, "lon"); lat <- ncvar_get(nc, "lat")
nx  <- length(lon); ny_lat <- length(lat)
yrs <- 1901:2021; ny <- length(yrs)
stopifnot(nc$dim[["time"]]$len == ny * 12)

annual <- array(NA_real_, dim = c(nx, ny_lat, ny))
for (k in seq_len(ny)) {
  m <- ncvar_get(nc, "triver",
                 start = c(1, 1, (k - 1) * 12 + 1),
                 count = c(nx, ny_lat, 12))          # [lon,lat,12]
  am <- rowMeans(m, dims = 2, na.rm = TRUE)          # mean over month dim
  am[is.nan(am)] <- NA
  annual[, , k] <- am
}
nc_close(nc)

lon_d <- ncdim_def("lon", "degrees_east", lon)
lat_d <- ncdim_def("lat", "degrees_north", lat)
tim_d <- ncdim_def("time", "years", yrs)
# compression=5 -> NetCDF4 deflate, shrinks this sparse field from ~120MB to ~22MB
var_d <- ncvar_def("rivtemp", "K", list(lon_d, lat_d, tim_d), missval = NA,
                   prec = "float", compression = 5)
fout  <- "input/rivertemp/watergap2-2e_20crv3-era5_obsclim_histsoc_default_rivtemp_global_annual_1901_2021.nc"
nco <- nc_create(fout, var_d, force_v4 = TRUE); ncvar_put(nco, var_d, annual); nc_close(nco)
cat("wrote", fout, "\n dims:", paste(dim(annual), collapse = "x"),
    "| valid river cells:", sum(!is.na(annual[, , 60])), "\n")

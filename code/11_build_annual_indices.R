# 11_build_annual_indices.R  -- SUPPLEMENTARY: annual-mean (calendar-year) oscillation indices,
# the year-round counterpart to the season-specific indices used in the main analysis.
# ENSO = mean of all 3-month ONI seasons in the year; NAO, IOD = mean of Jan..Dec.
# Output: input/clim_index/annual_indices.csv, with the SAME column names as
# seasonal_indices.csv (enso_djf, nao_djfm, iod_son) so the footprint engine (03/04) can be
# re-run on it unchanged via FP_INDEX_FILE (values are annual means, not seasonal).
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
idx_dir <- "input/clim_index/"; year_range <- 1950:2021

## ENSO: annual mean of the ONI (all overlapping 3-month seasons within each year)
oni <- read.table(paste0(idx_dir, "ONI"), header = TRUE)
enso <- aggregate(ANOM ~ YR, data = oni, FUN = mean, na.rm = TRUE)
names(enso) <- c("year", "enso_djf")

## NAO: annual mean (Jan..Dec)
nao_raw <- read.table(paste0(idx_dir, "NAO")); names(nao_raw) <- c("year", month.abb)
nao <- data.frame(year = nao_raw$year, nao_djfm = rowMeans(nao_raw[, month.abb], na.rm = TRUE))

## IOD: annual mean (Jan..Dec)
iod_raw <- read.table(paste0(idx_dir, "IOD")); names(iod_raw) <- c("year", month.abb)
iod <- data.frame(year = iod_raw$year, iod_son = rowMeans(iod_raw[, month.abb], na.rm = TRUE))

seasonal <- Reduce(function(a, b) merge(a, b, by = "year", all = TRUE), list(enso, nao, iod))
seasonal <- seasonal[seasonal$year %in% year_range, ]
seasonal <- seasonal[order(seasonal$year), ]
for (v in c("enso_djf", "nao_djfm", "iod_son")) seasonal[[v]] <- as.numeric(scale(seasonal[[v]]))

out <- paste0(idx_dir, "annual_indices.csv")
write.csv(seasonal, out, row.names = FALSE, quote = FALSE)
cat("Wrote", out, "-", nrow(seasonal), "years (", min(seasonal$year), "-", max(seasonal$year), ")\n")
cat("Correlation of annual vs season-specific indices (1950-2021):\n")
ss <- read.csv(paste0(idx_dir, "seasonal_indices.csv")); ss <- ss[match(seasonal$year, ss$year), ]
for (v in c("enso_djf","nao_djfm","iod_son"))
  cat(sprintf("  %-9s r = %+.2f\n", v, cor(seasonal[[v]], ss[[v]], use = "complete.obs")))

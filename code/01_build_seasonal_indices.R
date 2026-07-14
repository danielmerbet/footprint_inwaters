# build_seasonal_indices.R
# Build season-specific oscillation indices for the variance-partitioning analysis.
#   ENSO -> DJF   (ONI 3-month season, boreal winter peak)
#   NAO  -> DJFM  (Dec[y-1]+Jan+Feb+Mar, winter mode)
#   IOD  -> SON   (Sep+Oct+Nov, autumn Indian Ocean Dipole peak)
# Output: input/clim_index/seasonal_indices.csv  (year, enso_djf, nao_djfm, iod_son), 1950-2021.

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
idx_dir <- "input/clim_index/"
year_range <- 1950:2021

## ---- ENSO: ONI, DJF season ----------------------------------------------
# File columns: SEAS YR TOTAL ANOM  (ANOM = ONI value). DJF[y] = Dec[y-1]+Jan[y]+Feb[y].
oni <- read.table(paste0(idx_dir, "ONI"), header = TRUE)
enso <- oni[oni$SEAS == "DJF", c("YR", "ANOM")]
names(enso) <- c("year", "enso_djf")

## ---- NAO: monthly table, DJFM winter mean -------------------------------
# Columns: year Jan Feb ... Dec (13 cols). DJFM[y] = mean(Dec[y-1], Jan[y], Feb[y], Mar[y]).
nao_raw <- read.table(paste0(idx_dir, "NAO"))
names(nao_raw) <- c("year", month.abb)                 # year, Jan..Dec
nao <- data.frame(year = nao_raw$year)
dec_prev <- nao_raw$Dec[match(nao$year - 1, nao_raw$year)]   # December of previous year
nao$nao_djfm <- rowMeans(cbind(dec_prev, nao_raw$Jan, nao_raw$Feb, nao_raw$Mar), na.rm = TRUE)

## ---- IOD: monthly table, SON autumn mean --------------------------------
# Columns: year Jan..Dec (13 cols). SON[y] = mean(Sep[y], Oct[y], Nov[y]).
iod_raw <- read.table(paste0(idx_dir, "IOD"))
names(iod_raw) <- c("year", month.abb)
iod <- data.frame(year = iod_raw$year,
                  iod_son = rowMeans(iod_raw[, c("Sep", "Oct", "Nov")], na.rm = TRUE))

## ---- Merge, subset, standardize -----------------------------------------
seasonal <- Reduce(function(a, b) merge(a, b, by = "year", all = TRUE),
                  list(enso, nao, iod))
seasonal <- seasonal[seasonal$year %in% year_range, ]
seasonal <- seasonal[order(seasonal$year), ]

# z-score each index (comparable regression coefficients; does not affect R^2)
for (v in c("enso_djf", "nao_djfm", "iod_son")) {
  seasonal[[v]] <- as.numeric(scale(seasonal[[v]]))
}

out <- paste0(idx_dir, "seasonal_indices.csv")
write.csv(seasonal, out, row.names = FALSE, quote = FALSE)

cat("Wrote", out, "-", nrow(seasonal), "years (",
    min(seasonal$year), "-", max(seasonal$year), ")\n")
cat("NA counts:\n"); print(colSums(is.na(seasonal)))
cat("Correlations among indices (collinearity check):\n")
print(round(cor(seasonal[, -1], use = "pairwise.complete.obs"), 2))
cat("Summary:\n"); print(summary(seasonal[, -1]))

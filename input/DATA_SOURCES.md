# Data sources & provenance

Full provenance for every input to *The footprint of climate oscillations in global inland waters*.
All gridded data are on a common **0.5° regular lon/lat grid, CRS EPSG:4326**. Water-temperature
and meteorological variables are annual means over **1950–2021** (built from 1901–2021 source files).

> **Access model.** The large NetCDF inputs below are third-party ISIMIP model output and reanalysis
> with their own repositories, DOIs and terms of use — download them from the sources listed rather
> than from this code repository. The oscillation indices are small text files and are included.

---

## 1. Water temperature (model output, ISIMIP3a)

| Property | Lake temperature | River temperature |
|---|---|---|
| Variables | `surftemp`, `bottemp` (lake surface / bottom water temperature, K) | `triver` → aggregated to annual `rivtemp` (river water temperature, K) |
| Model | GOTM (1-D lake model), ISIMIP Lake Sector | WaterGAP2-2e global hydrological model |
| Reference | Golub et al. (2022), *Geosci. Model Dev.* 15, 4597–4623, doi:10.5194/gmd-15-4597-2022 | Müller Schmied et al. (2021), *Geosci. Model Dev.* 14, 1037–1079, doi:10.5194/gmd-14-1037-2021 |
| Forcing | 20CRv3-ERA5 (also 20CRv3-W5E5, GSWP3-W5E5 for multi-forcing test, `code/09`) | 20CRv3-ERA5 |
| Files | `laketemp/gotm_<forcing>_obsclim_histsoc_default_{surftemp,bottemp}_global_annual_*.nc` | `rivertemp/watergap2-2e_..._{triver,rivtemp}_global_annual_1901_2021.nc` |
| Format / units | NetCDF-4 (CF), kelvin | NetCDF-4 (CF), kelvin |
| Access | ISIMIP Repository — <https://data.isimip.org/> (ISIMIP3a, Lake / Water Global sectors) | ISIMIP Repository — <https://data.isimip.org/> (ISIMIP3a, Water Global sector) |

*Note:* the monthly `triver` source is aggregated to annual `rivtemp` by `code/02_build_annual_rivertemp.R`.

---

## 2. Meteorological forcing (local climate predictors)

| Property | Value |
|---|---|
| Variables | `tas` (near-surface air temperature), `rsds` (downwelling shortwave), `rlds` (downwelling longwave), `sfcwind` (surface wind speed), `pr` (precipitation) — plus `tasmax`, `tasmin`, `hurs`, `ps` (available, not used) |
| Dataset | 20CRv3-ERA5, ISIMIP3a obsclim protocol |
| Reference | Frieler et al. (2024), *Geosci. Model Dev.* 17, 1–51, doi:10.5194/gmd-17-1-2024 |
| Files | `clim_data/20crv3-era5_obsclim_<var>_global_yearly_1901_2021.nc` |
| Format | NetCDF-4 (CF); SI units per variable |
| Access | ISIMIP Repository — <https://data.isimip.org/> (ISIMIP3a input, climate) |

---

## 3. Oscillation indices (included in `clim_index/`)

Season-specific indices are built by `code/01_build_seasonal_indices.R` into
`clim_index/seasonal_indices.csv` (columns `year, enso_djf, nao_djfm, iod_son`, standardised).
Raw source files and their URLs (`clim_index/reference`):

| Index | Definition used | Provider | Source URL |
|---|---|---|---|
| ENSO | Oceanic Niño Index (ONI), DJF | NOAA Climate Prediction Center | <https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt> |
| NAO | NAO index, DJFM | NOAA Climate Prediction Center | <https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii.table> |
| IOD | Dipole Mode Index (DMI, HadISST-based), SON | NOAA Physical Sciences Laboratory | <https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/dmi.had.long.data> |

*(The `clim_index/` folder also contains SOI, PDO, NPO, NPGO for reference; only ONI, NAO and DMI
enter the analysis. ENSO uses ONI, not SOI.)*

---

## 4. Lake bathymetry

| Property | Value |
|---|---|
| Variables | `max_lakedepth` (maximum lake depth, m), `surface_area` (m²) |
| Files | `lake_bathymetry/max_lakedepth.nc`, `surface_area.nc` |
| Use | `max_lakedepth` used by `code/08_depth_memory.R` (direct-fraction vs depth) |
| Source | ISIMIP Lake Sector morphometry (Golub et al., 2022) — <https://data.isimip.org/> |
| Format | NetCDF, 0.5°, EPSG:4326 |

---

## 5. Ancillary (present but not used in the current analysis)

| File | Note |
|---|---|
| `clim_data/koppen_geiger_0p00833333_1991-2020.tif` | Köppen–Geiger climate zones (Beck et al., 2018). Retained in the input folder but **not used** by the committed analysis. |

---

## Licences / terms of use
Input datasets retain the licences of their original providers. ISIMIP data are distributed under
the ISIMIP terms of use (<https://www.isimip.org/gettingstarted/terms-of-use/>); NOAA indices are
in the public domain. Cite the references above when reusing. The **derived** products in `output/`
are released under CC BY 4.0 (see `../LICENSE-CC-BY-4.0.md`).

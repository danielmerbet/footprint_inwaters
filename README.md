# The footprint of climate oscillations in global inland waters

[![License: MIT](https://img.shields.io/badge/Code-MIT-blue.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/Data-CC%20BY%204.0-lightgrey.svg)](LICENSE-CC-BY-4.0.md)
[![DOI](https://img.shields.io/badge/Zenodo-DOI%20pending-orange.svg)](https://doi.org/FILL-ZENODO-DOI)

Reproducible analysis quantifying and explaining how the leading modes of interannual
climate variability — **ENSO, IOD and NAO** — explain year-to-year temperature variability in
the world's **rivers** and **lakes (surface & bottom)**, on a common 0.5° global grid, 1950–2021.

The headline metric is the oscillation **footprint** *F* — the fraction of detrended interannual
water-temperature variance jointly explained by the three modes,
*F* = adjusted *R²*(*T′* ~ ENSO + IOD + NAO). Grid cells where *F* is field-significant
(Benjamini–Hochberg FDR < 0.05) are termed **teleconnected**.

> **Manuscript:** *The footprint of climate oscillations in global inland waters* —
> D. Mercado-Bettín & R. Marcé (CEAB-CSIC), in preparation for submission.

---

## Contents
- [Highlights](#highlights)
- [Repository structure](#repository-structure)
- [Software environment](#software-environment)
- [Reproduce the analysis](#reproduce-the-analysis)
- [Data & provenance](#data--provenance)
- [Outputs](#outputs)
- [Interactive artifacts](#interactive-artifacts)
- [How to cite](#how-to-cite)
- [License](#license)
- [FAIR + Reproducibility statement](#fair--reproducibility-statement)

---

## Highlights
- **Teleconnected extent / median *F*:** river 3.0% / 0.28 · lake surface 7.4% / 0.30 · lake bottom 8.7% / 0.25.
- **Dominant mode:** ENSO everywhere (72% river, 78% lake surface); NAO rises to 26% at the lake bottom.
- **All compartments respond at lag 0** (contemporaneous; tropical-ENSO composite *r* ≈ 0.82–0.85).
- **Mechanism:** > 97% atmosphere-mediated at surfaces; only the deep lake develops a direct /
  memory pathway (~20% of teleconnected bottom cells), which scales with lake depth
  (Spearman *ρ* = 0.38, *p* < 10⁻¹⁰⁰).
- **Field significance & robustness:** 500-surrogate phase-randomization null (*p* < 0.002 for
  extent and strength); footprint replicates across three meteorological forcings (pixelwise
  *r* = 0.77–0.88); block-CV *R²* 0.21 (rivers) / 0.23 (lakes) with ~100% positive out-of-sample
  skill; index VIF ≈ 1; *F* insensitive to lag and to detrend/first-difference band isolation.

---

## Repository structure
```
footprint_inwaters/
├── README.md                     ← this file
├── CITATION.cff                  ← how to cite (GitHub + Zenodo)
├── .zenodo.json                  ← Zenodo deposit metadata
├── codemeta.json                 ← CodeMeta software metadata (FAIR)
├── LICENSE                       ← MIT (code)
├── LICENSE-CC-BY-4.0.md          ← CC BY 4.0 (data, figures, text)
├── install.R                     ← R dependencies & environment
│
├── input/
│   ├── DATA_SOURCES.md           ← full data provenance table (sources, DOIs, licences)
│   ├── clim_index/               ← oscillation indices: seasonal_indices.csv (+ raw ONI/NAO/IOD, reference)
│   ├── clim_data/                ← 20CRv3-ERA5 local climate (tas, rsds, rlds, sfcwind, pr, …) [third-party]
│   ├── laketemp/                 ← GOTM lake surface & bottom temperature, ISIMIP3a (3 forcings)
│   ├── rivertemp/                ← WaterGAP2-2e river temperature (annual rivtemp & triver)
│   └── lake_bathymetry/          ← max_lakedepth.nc, surface_area.nc
│
├── code/
│   ├── 01_build_seasonal_indices.R   raw indices → seasonal_indices.csv
│   ├── 02_build_annual_rivertemp.R   monthly triver → annual rivtemp (memory-safe, deflate)
│   ├── 03_footprint.R                footprint F, per-index R², cluster lags, eff-DOF sig, BH-FDR
│   ├── 04_mechanism.R                two-block variation partitioning: direct vs mediated + mediator
│   ├── 05_robustness.R               block cross-validation, VIF, lag/detrend sensitivity
│   ├── 06a_export_interactive_lakes.R  } compact JSON for the interactive explorers
│   ├── 06b_export_interactive_rivers.R }
│   ├── 07_surrogate_null.R           phase-randomization field-significance + selection null
│   ├── 08_depth_memory.R             direct/memory fraction vs maximum lake depth
│   ├── 09_multiforcing.R             footprint replication across 3 meteorological forcings
│   ├── 10_supp_stats.R               index collinearity, negative VP components, first-difference
│   ├── 11_build_annual_indices.R     calendar-year annual-mean indices (season-specific sensitivity)
│   └── figures/
│       ├── F1_phenomenon.R  F2_footprint.R  F3_mechanism.R  F4_synthesis.R   (main figures)
│       └── FS_nao_iod.R      NAO & IOD analogues of F1 and F4b/d (supplementary)
│
└── output/
    ├── footprint/{rivtemp,surftemp,bottemp}/          GeoTIFFs: F, p, q, per-index R²/lag/sig
    │   └── …/mechanism/                               direct/mediated/local R², mediated_frac, mediator
    └── figures/
        ├── F1_phenomenon.png … F4_synthesis.png       main manuscript figures
        ├── interactive/                               self-contained HTML explorers (see below)
        └── supplementary/                             robustness/, NAO-IOD figures, legacy figures
```
*Excluded from version control (`.gitignore`): `largefiles/` (uncompressed source NetCDFs) and
`manuscript/` (paper drafts, released separately on acceptance).*

---

## Software environment
- **R ≥ 4.1** (developed on R 4.1.2). Analysis packages:
  `terra` (≥1.8), `ncdf4` (≥1.23), `sf` (≥1.0), `ggplot2` (≥3.5), `patchwork`, `jsonlite`,
  `tidyterra`, `rnaturalearth`, `viridis`, `ggrepel`, `ggnewscale`.
- **System libraries** (for `terra`/`sf`/`ncdf4`): GDAL, GEOS, PROJ, and the NetCDF C library.
- **Optional:** `cdo` (Climate Data Operators) for the NetCDF compression used in data preparation;
  a TeX distribution with `pdflatex` to build the manuscript.

Install the R dependencies with:
```r
source("install.R")
```
The exact tested versions are recorded in [`install.R`](install.R). All scripts locate the
repository root automatically (they walk up to the folder containing `input/` and `code/`), so
they run unchanged via `Rscript`, from any working directory, or sourced in RStudio.

---

## Reproduce the analysis
All paths are relative to the repository root; no editing required. Recommended order:

| Step | Command | Produces |
|---|---|---|
| Build indices | `Rscript code/01_build_seasonal_indices.R` | `input/clim_index/seasonal_indices.csv` |
| Build river temp | `Rscript code/02_build_annual_rivertemp.R` | annual `rivtemp` NetCDF |
| **Footprint** | `Rscript code/03_footprint.R` | `output/footprint/*/` GeoTIFFs |
| **Mechanism** | `Rscript code/04_mechanism.R` | `output/footprint/*/mechanism/` |
| Robustness | `Rscript code/05_robustness.R` | `output/figures/supplementary/robustness/` |
| Surrogate null | `Rscript code/07_surrogate_null.R` | field-significance null (Fig S2) |
| Depth–memory | `Rscript code/08_depth_memory.R` | depth vs direct-fraction (Fig S1) |
| Multi-forcing | `Rscript code/09_multiforcing.R` | forcing replication (Fig S3) |
| Supp. statistics | `Rscript code/10_supp_stats.R` | index collinearity, VP components, first-difference |
| Explorer data | `Rscript code/06a_export_interactive_lakes.R` · `…/06b_…rivers.R` | explorer JSON |
| Main figures | `Rscript code/figures/F1_phenomenon.R` … `F4_synthesis.R` | `output/figures/F1–F4.png` |
| NAO/IOD figures | `Rscript code/figures/FS_nao_iod.R` | `output/figures/supplementary/FS_*.png` |

`03` → `04` → `05`/`07`–`10` is the core chain; the figure scripts read the products of `03`/`04`.
Steps `01`/`02` rebuild inputs that are already provided.

**Sensitivity variant (annual-mean indices).** `03` and `04` accept two environment variables —
`FP_INDEX_FILE` (index CSV) and `FP_OUT_ROOT` (output tree; defaults reproduce the main analysis).
To reproduce the annual-mean-index supplement (Figs S9–S10):
```sh
Rscript code/11_build_annual_indices.R
FP_INDEX_FILE=input/clim_index/annual_indices.csv FP_OUT_ROOT=output/footprint_annual/ Rscript code/03_footprint.R
FP_INDEX_FILE=input/clim_index/annual_indices.csv FP_OUT_ROOT=output/footprint_annual/ Rscript code/04_mechanism.R
FP_OUT_ROOT=output/footprint_annual/ FIG_OUT=output/figures/supplementary/FS_annual_footprint.png Rscript code/figures/F2_footprint.R
FP_OUT_ROOT=output/footprint_annual/ FIG_OUT=output/figures/supplementary/FS_annual_mechanism.png Rscript code/figures/F3_mechanism.R
```

**Method in brief.** Annual water-temperature and index series are linearly detrended to isolate
the interannual band; season-specific indices (ENSO = DJF, NAO = DJFM, IOD = SON) are standardised
(mean 0, SD 1). Per cell, *F* = adjusted *R²* of the multiple regression of *T′* on the three
(cluster-lagged) indices; significance uses an autocorrelation-corrected (effective degrees of
freedom) *F*-test, with Benjamini–Hochberg FDR for field significance. Mechanism uses two-block
adjusted-*R²* variation partitioning (Borcard/Legendre) into direct vs locally-mediated pathways.
An interactive walkthrough of the *F* computation is provided (see [Interactive artifacts](#interactive-artifacts)).

---

## Data & provenance
Full machine- and human-readable provenance — sources, versions, DOIs, formats, spatial/temporal
resolution and licences — is in **[`input/DATA_SOURCES.md`](input/DATA_SOURCES.md)**. Summary:

| Dataset | Variable(s) | Source | Format | Grid / period |
|---|---|---|---|---|
| Lake temperature | surface, bottom | GOTM, ISIMIP Lake Sector (Golub et al., 2022) | NetCDF | 0.5°, 1901–2021 |
| River temperature | `triver` → `rivtemp` | WaterGAP2-2e (Müller Schmied et al., 2021) | NetCDF | 0.5°, 1901–2021 |
| Meteorological forcing | tas, rsds, rlds, sfcwind, pr | 20CRv3-ERA5, ISIMIP3a (Frieler et al., 2024) | NetCDF | 0.5°, 1901–2021 |
| ENSO index | ONI (DJF) | NOAA CPC | ASCII/CSV | monthly |
| NAO index | DJFM | NOAA CPC | ASCII/CSV | monthly |
| IOD index | DMI (SON) | NOAA PSL (HadISST-based) | ASCII/CSV | monthly |
| Lake bathymetry | max depth, area | ISIMIP Lake Sector | NetCDF | 0.5° |

**Two-tier FAIR deposit.** Raw meteorological and temperature NetCDFs (~0.8 GB) are third-party
ISIMIP/reanalysis products with their own access terms and identifiers; obtain them from the
[ISIMIP Repository](https://data.isimip.org/) and the NOAA URLs in `input/DATA_SOURCES.md`.
This code repository (GitHub) contains the analysis code, small derived products and metadata; the
**derived data products and figures are archived with a persistent DOI on Zenodo** (see below).

---

## Outputs
- `output/footprint/{rivtemp,surftemp,bottemp}/` — per-cell GeoTIFFs: footprint `F`
  (`footprint_adjR2.tiff`), effective-DOF `p` and BH `q`, per-index adjusted `R²`, lag and
  significance masks, and (under `mechanism/`) direct/mediated/local `R²`, mediated fraction and
  dominant mediator.
- `output/figures/F1–F4.png` — the four-figure argument:
  **F1 Phenomenon** (ENSO± composites) → **F2 Footprint** (*F* + dominant index) →
  **F3 Mechanism** (delivering pathway) → **F4 Synthesis** (exposed vs buffered).
- `output/figures/supplementary/` — robustness panels, NAO/IOD analogues (Figs S5–S8), legacy figures.

---

## Interactive artifacts
Self-contained, dependency-free HTML (open in any browser):
- `output/figures/interactive/rivers-explorer.html`, `lakes-explorer.html` — per-cell explorers of
  the footprint, indices and model fit.
- `output/figures/interactive/footprint_method_explainer.html` — a step-by-step, interactive
  walkthrough that recomputes *F* live on real example cells (detrend, standardise, regress,
  adjusted *R²*, significance), reproducing the pipeline values exactly.

---

## How to cite
If you use this code, data or method, please cite both the software/data archive and (when
published) the article. Citation metadata is in [`CITATION.cff`](CITATION.cff); GitHub renders a
"Cite this repository" button from it and Zenodo mints a DOI on release.

```
Mercado-Bettín, D., & Marcé, R. (2026). The footprint of climate oscillations in global inland
waters: analysis code and derived data (v1.0.0) [Software/Dataset]. Zenodo.
https://doi.org/FILL-ZENODO-DOI
```

---

## License
- **Code** (`code/`, `install.R`): [MIT](LICENSE).
- **Data products, figures and text** (`output/`, derived data, manuscript): [CC BY 4.0](LICENSE-CC-BY-4.0.md).
- Third-party input datasets retain the licences of their original providers (see `input/DATA_SOURCES.md`).

---

## FAIR + Reproducibility statement
This repository is designed to meet the **FAIR** principles and full computational reproducibility:

- **Findable** — rich metadata (`CITATION.cff`, `.zenodo.json`, `codemeta.json`), descriptive
  keywords, and a persistent **Zenodo DOI** minted from a tagged GitHub release.
- **Accessible** — open, free formats retrievable over HTTP(S); explicit access routes for every
  dataset (`input/DATA_SOURCES.md`); no proprietary software required.
- **Interoperable** — standard, self-describing formats (CF-NetCDF, GeoTIFF, CSV, JSON, HTML) and
  documented variables, units and CRS (EPSG:4326).
- **Reusable** — clear licences (MIT + CC BY 4.0), complete provenance, and pinned software
  versions.
- **Reproducible (+R)** — the full analysis regenerates from raw inputs with the numbered scripts
  above; all code is path-agnostic (auto-detects the repository root) and runs on open-source R.

**Contact:** Daniel Mercado-Bettín — `daniel.mercado@csic.es` · Centre d'Estudis Avançats de Blanes
(CEAB-CSIC), Blanes, Spain.

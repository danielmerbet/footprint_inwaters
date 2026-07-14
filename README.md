# The footprint of climate oscillations in global inland waters

Quantifying and explaining how the leading modes of interannual climate variability —
**ENSO, IOD and NAO** — organise year-to-year temperature variability in the world's **lakes
(surface & bottom)** and **rivers**, on a common 0.5° global grid, 1950–2021.

The headline metric is the oscillation **footprint** *F* = the fraction of detrended interannual
water-temperature variance jointly explained by the three modes. Cells where *F* is
field-significant (BH-FDR < 0.05) are termed **teleconnected**.

## Folder layout
```
footprint_inwaters/
  input/
    clim_index/     season-specific indices (seasonal_indices.csv) + raw sources
    clim_data/      20CRv3-ERA5 local climate (tas,rsds,rlds,sfcwind,pr) + Köppen-Geiger
    laketemp/       GOTM lake surface & bottom temperature (ISIMIP3a)
    rivertemp/      WaterGAP2-2e river temperature (triver monthly + annual rivtemp)
    lake_bathymetry/
  code/
    01_build_seasonal_indices.R     raw indices -> seasonal_indices.csv
    02_build_annual_rivertemp.R     monthly triver -> annual rivtemp (memory-safe)
    03_footprint.R                  footprint F, per-index R², cluster lags, FDR  (all 3 compartments)
    04_mechanism.R                  two-block partition: direct vs mediated + dominant mediator
    05_robustness.R                 block-CV, VIF, lag/detrend sensitivity
    06a_export_interactive_lakes.R  } JSON for the interactive explorers
    06b_export_interactive_rivers.R }
    figures/  F1_phenomenon.R  F2_footprint.R  F3_mechanism.R  F4_synthesis.R
  output/
    footprint/{surftemp,bottemp,rivtemp}/   GeoTIFFs (F, p, q, per-index R²/lag/sig, mechanism/)
    figures/  F1..F4 .png
              interactive/   lakes-explorer.html, rivers-explorer.html (self-contained)
              supplementary/ robustness/, S1_phase_composites_atlas.png, S2_catboost_importance.png
  manuscript/
    footprint_inwaters_GRL_draft.md   full GRL draft (target: Geophysical Research Letters)
```

## Run order
`03 → 04 → 05` regenerate all footprint/mechanism/robustness products; `figures/F1..F4` build the
manuscript figures; `06a/06b` regenerate the explorer data. `01`/`02` rebuild inputs (already
provided). Every script sets `base <- ".../footprint_inwaters"` and uses paths relative to it.

## The four-figure argument
**F1 Phenomenon** (ENSO± composites, sign-flip) → **F2 Footprint** (F + dominant index, three
compartments) → **F3 Mechanism** (which local pathway delivers the signal) → **F4 Synthesis**
(exposed surface waters vs the buffered deep lake).

## Headline results
- **Teleconnected extent / median F:** lake surface 7.4% / 0.30 · lake bottom 8.7% / 0.25 · river 3.0% / 0.28.
- **Dominant mode:** ENSO everywhere (lake-surface 78%, river 72%); NAO rises to 26% at the lake bottom.
- **All compartments respond at lag 0** (contemporaneous; tropical-ENSO composite r ≈ 0.82–0.85).
- **Mechanism:** >97% atmosphere-mediated at surfaces; only the deep lake develops a direct/memory
  pathway (~20% of teleconnected bottom cells). Precipitation is a river-specific mediator (peaks
  13% in temperate rivers).
- **Robustness:** block-CV R² 0.23 (lakes) / 0.21 (rivers), ~100% positive out-of-sample skill;
  index VIF ≈ 1; F insensitive to lag & detrend order.

## Data provenance
Lake (GOTM) and river (WaterGAP2-2e) temperatures are ISIMIP3a simulations forced by the
20CRv3-ERA5 reanalysis. Climate indices: NOAA CPC (ENSO/ONI, NAO) and NOAA PSL (IOD/DMI).

## Supplementary (legacy) figures
`S1_phase_composites_atlas.png` (phase-composite anomaly maps) and `S2_catboost_importance.png`
(CatBoost predictor importance) come from an earlier analysis generation, retained as supplement;
their generating code was retired in favour of the detrended variance-partitioning pipeline above.

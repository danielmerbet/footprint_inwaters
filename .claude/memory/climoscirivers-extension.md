---
name: climoscirivers-extension
description: Rivers extension of the aquatic-sentinels analysis (ClimOsciRivers repo)
metadata: 
  node_type: memory
  type: project
  originSessionId: d5951ec6-1e37-4e52-8b0d-218beb2d185f
---

**MERGED (2026-07-13):** rivers are no longer a separate repo — everything consolidated into `/home/dmercado/Documents/footprint_inwaters/` (see [[climoscilakes-paper-plan]]). River temp is the `rivtemp` compartment: input at `input/rivertemp/`, outputs at `output/footprint/rivtemp/`, built by the same unified `code/03_footprint.R` / `04_mechanism.R`. The old `ClimOsciRivers/` folder is DELETED. Paths below are historical.

Extension of [[climoscilakes-paper-plan]] from lakes to **rivers**, broadening the framing to inland waters. (Originally separate repo `/home/dmercado/Documents/ClimOsciRivers/`.)

**Data:** river temperature = WaterGAP2-2e `triver` (ISIMIP3a, 20crv3-era5 obsclim, monthly 1901–2021, 269MB). NOTE: I initially wrongly said river temp didn't exist; user was right — it's `triver` in the ISIMIP water-global sector. `build_annual_rivtemp.R` aggregates monthly→annual `rivtemp` NetCDF, **year-by-year** via `ncvar_get(start/count)` to avoid OOM (loading all 1452 months at once was killed, exit 137). `input/clim_index` + `input/clim_data` are symlinks to the lakes repo (shared indices + meteorology).

**Pipeline mirrors lakes** (same method, grid, 1950–2021): `variance_partitioning.R` → `fig3_river_sentinel_map.R` (single panel) → `fig3b_river_dominant_index.R` → `mediation_partitioning.R` → `fig4_river_mechanism.R` → `robustness_checks.R`. Mechanism local block **adds precipitation** (`locvars` = tas,rsds,rlds,sfcwind,**pr**) — river-specific hydrological pathway; fig4 has 6 categories.

**Results:** 63,978 river cells, 3.0% FDR-significant, median F=0.28. Dominant: ENSO 72/IOD 16/NAO 12%. Footprint on Amazon, Maritime Continent, E Africa, Iberia, Pacific NW. Median mediated fraction **0.98** (air-temp 61/longwave 18/wind 9/solar 7/precip 6%; precip peaks 13% in temperate rivers). Robustness: CV R² 0.208 vs 0.280 in-sample, 99.8% positive skill; VIF fine.

**Comparative headline (rivers vs lakes):** same per-site strength (F≈0.28) but rivers are **fast integrators** — far fewer teleconnected (3% vs 7–9%) and ~no thermal memory (98% mediated), vs lakes (esp. deep) as **slow integrators** with a growing direct/memory pathway. Full write-up in repo `README.md`.

**Interactive river explorer** built (`export_interactive_data.R` → 2.3MB JSON → `output/figures/fig3/river-sentinels-explorer.html`, single-layer adaptation of the lakes explorer) and **published as Artifact**: https://claude.ai/code/artifact/61d2bc3d-9975-4da3-a8fc-4b1d58461701 .

**Combined "aquatic systems" manuscript** drafted: `ClimOsciLakes/manuscript/ClimOsciAquatic_GRL_draft.md` (lakes-only draft kept alongside). Full rivers pipeline complete (analysis+figs+interactive+robustness), mirroring lakes.

**Figure arc restructured (2026-07-13):** user found the old 4 figures overlapping/unconvincing. Diagnosis: old Fig 1 (phase-composite anomaly maps, `output/figures/fig1/`) and Fig 2 (CatBoost variable-importance bars, `output/figures/fig2/`) are the OLD-method predecessors of the new Fig 3 (footprint) & Fig 4 (mechanism) — two generations of the same two questions. The user's own README `notas` already condemned the CatBoost approach. New arc built in `manuscript/figs_aquatic/` (scripts F1–F4, PNGs in `out/`): **F1 phenomenon** (detrended ENSO± composites, sign-flip, lakes+rivers) → **F2 footprint** (F + dominant index; lake-surf/bot/river) → **F3 mechanism** (replaces CatBoost) → **F4 synthesis** (strength-vs-extent, memory axis, mediator identity, lag-0 composite). Old Fig 1/2 demoted to supplement. Framing corrected to **exposed vs buffered** (NOT fast/slow — lag is 0 everywhere; see [[climoscilakes-paper-plan]] correction).

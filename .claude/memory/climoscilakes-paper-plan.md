---
name: climoscilakes-paper-plan
description: Locked scope & framing for the ClimOsciLakes lake-sentinel paper
metadata: 
  node_type: memory
  type: project
  originSessionId: d5951ec6-1e37-4e52-8b0d-218beb2d185f
---

**CONSOLIDATED (2026-07-13):** the whole project now lives in ONE folder `/home/dmercado/Documents/footprint_inwaters/` (old `ClimOsciLakes` + `ClimOsciRivers` folders DELETED). Layout: `input/`, `code/` (01–06 numbered pipeline + `figures/F1–F4`), `output/footprint/{surftemp,bottemp,rivtemp}/`, `output/figures/`, `manuscript/footprint_inwaters_GRL_draft.md`. Manuscript title = **"The footprint of climate oscillations in global inland waters"**. Operational term is **teleconnected** cells (not "sentinel" — that word appears only in the intro nod to Adrian 2009/Williamson 2009). Figure arc = F1 phenomenon → F2 footprint → F3 mechanism → F4 synthesis. Manual QGIS polygons preserved under `output/figures/supplementary/legacy_qgis_clusters/`.

Building a high-impact paper: inland waters (lakes surface+bottom, rivers) and their interannual-oscillation footprint.

**Locked decisions (as of 2026-07-10/11):**
- Scope: lakes only, surface **and** bottom temperature (GOTM 20crv3-era5, 1950–2021).
- Oscillations: ENSO, IOD, NAO (interannual); PDO dropped (interdecadal).
- Headline framing: **commonality-based variance partitioning** — fraction of interannual lake thermal variance attributable to oscillations (Unique_osc + Common), mapped globally as a "sentinel-strength" map, surface vs bottom.
- CatBoost work repurposed as the mechanism/mediation layer, not the headline.
- Target journal: **GRL** (see [[journal-preference-nonprofit]]).

**Method decisions (reviewed 2026-07-11):**
- Headline footprint F simplifies algebraically to **adj-R²(lake temp ~ ENSO+IOD+NAO)** = R²_O (commonality's Unique_osc+Common). No local block or causal assumption needed for the headline; commonality (Borcard/Legendre variation partitioning, `vegan::varpart`) is only the Fig-4 mechanism layer (direct vs mediated), framed descriptively.
- **Detrend** lake temp to isolate interannual band before regressing (linear default, SSA sensitivity).
- **Season-specific indices** (ENSO≈DJF, IOD≈SON, NAO≈DJFM) vs annual lake-temp anomaly.
- **Cluster-informed fixed lag** (from aggregated CCF per cluster+index), not per-pixel selection.
- Baked-in defensibility: effective-DOF for autocorrelation, FDR field significance, index attribution via relative weights (ENSO/IOD collinear), lag choice inside the permutation null.

**Build progress (2026-07-11):**
- `build_seasonal_indices.R` + `variance_partitioning.R` built & validated (both depths). Footprint lands on correct teleconnection regions (ENSO→Peru/SE Asia; NAO→Europe/N Atlantic; IOD→Indian Ocean rim).
- **Key finding:** raw IOD "significance" was ~83% spurious shared warming trend (IOD r=+0.48 with year). Fix: the footprint pipeline is **self-contained & fully detrended** — detrends lake temp AND indices, computes its own per-index detrended significance, defines cluster-informed lag from coherent significant patches (≥5 px). Does NOT use the legacy raw correlation masks/clusters, so Fig 1/2 + hand-refined QGIS polygons stay untouched.
- FDR (BH) field significance in the engine. Headline: median F among FDR-significant lakes = **surf 0.30, bot 0.25**.
- **CORRECTION (2026-07-13):** the old "bottom lags 2–3 years / surface faster" claim is **FALSE** — the detrended engine assigns **lag 0 to ~98% of sentinels** (ENSO & NAO always 0; only IOD lag-1 in ~2%). All compartments (river, lake surface, lake bottom) respond **contemporaneously** (tropical-ENSO composite tracks ENSO at lag 0, r≈0.83 in all three). The surface-vs-deep contrast is **memory/amplitude, not timing**: median direct fraction surf 0.008 < river 0.024 < bottom 0.054. Framing is **exposed** (rivers+lake surfaces, fully mediated) **vs buffered** (deep lake stores memory). Do NOT reintroduce lag language.

**Figures built (all in `output/figures/`):** Fig 3 sentinel map (`fig3_sentinel_map.R`), Fig 3 companion dominant-index map (`fig3b_dominant_index.R`), interactive explorer (`lake-sentinels-explorer.html`, published Artifact), Fig 4 mechanism map (`mediation_partitioning.R`+`fig4_mechanism.R`).

**Fig 4 mechanism result:** two-block variation partitioning (adj-R² commonality, a+b=F). Surface footprint is 92% air-temperature-mediated; at depth the direct/lake-memory pathway grows to 21% (median mediated fraction 0.98→0.93; cold deep lakes 0.75); wind-mixing mediates in tropical maritime lakes.

**Robustness done** (`robustness_checks.R`): block-CV R² 0.227 vs in-sample 0.304 (100% lakes positive skill); index VIF 1.02–1.05; local VIF tas 4.3/rlds 5.3 (thermal-radiative collinearity), wind 1.2; median F insensitive to lag & detrend order. → `output/figures/supplementary/robustness/`.

**Manuscript draft done**: `manuscript/ClimOsciLakes_GRL_draft.md` — full GRL draft with live numbers, placeholder refs to fill in a reference manager.

**Status: full analysis + figures + interactive + robustness + draft complete.** Remaining is authorial: verify references, finalize author list/affiliations, polish prose, convert to journal template, and (optional) revisit Fig 1/2 detrending consistency if reviewers ask.

**LaTeX + revisions built (2026-07-13/14):** `manuscript/` now has `footprint_inwaters_GRL.tex/.pdf` (v1, 11pp, mathptmx/natbib/line-numbers, 4 figs embedded), `..._v2.tex/.pdf` (revised, 12pp), `..._SI.tex/.pdf` (Supporting Info, 5pp), and `..._trackchanges.tex/.pdf` (latexdiff v1→v2). **R. Marcé is 2nd author in all versions.** Rebuild: `pdflatex` twice per file; track-changes via `latexdiff footprint_inwaters_GRL.tex footprint_inwaters_GRL_v2.tex`. **v3 (2026-07-14):** `..._v3.tex/.pdf` + `..._v3_trackchanges.pdf` (latexdiff v2→v3). v3 adds **affiliation = Centre d'Estudis Avançats de Blanes (CEAB-CSIC), Blanes, Spain** (shared, single) and NAO/IOD supplementary figures + main-text discussion. Author order: Mercado-Bettín, Marcé, [co-authors].

**NAO & IOD supplementary figures (2026-07-14):** `code/figures/FS_nao_iod.R` → `output/figures/supplementary/FS_{nao,iod}_{phenomenon,synthesis}.png` (SI Figs S5–S8). Phenomenon = F1 analogue (2×2 sign-flipping composites, own index quartiles + `{index}_sig.tiff`); synthesis = F4b memory-axis + F4d contemporaneous response over that index's DOMINANT cells, **sign-aligned** (NAO/IOD have spatial dipoles so raw areal mean cancels — ENSO tropics doesn't). Results: both modes lag-0 (CCF peak 0; NAO r=0.64–0.78, IOD r=0.77–0.78); memory axis holds — deep-lake direct fraction NAO 0.31 (winter convective mixing reaches depth), IOD 0.09, vs ≤0.05 surface/river. Main text now argues exposed-vs-buffered is general, not ENSO-only.

**Referee-facing revisions implemented** (new scripts `code/07_surrogate_null.R`, `08_depth_memory.R`, `09_multiforcing.R`, `10_supp_stats.R`; outputs in `output/figures/supplementary/robustness/`): (1) **phase-randomized field-significance null** (500 surrogates): obs 3–9% teleconnected vs null ~0%, p<0.002; selection-aware strength obs medF 0.25–0.30 vs null top-N 0.07–0.10, p<0.002. (2) **depth–memory**: direct fraction rises with max lake depth at BOTTOM (Spearman +0.38, p<1e-100; 0.02→0.34 over 3–100m) but flat at surface (p=0.46) — uses the previously-unused `input/lake_bathymetry/max_lakedepth.nc`. (3) **multi-forcing** (20crv3-era5/20crv3-w5e5/gswp3-w5e5, 1950–2019): median F stable 0.24–0.30, pixel r(F) 0.77–0.88 vs ERA5, Jaccard 0.25–0.54 (strength/pattern robust, extent forcing-sensitive). (4) index cross-corr |r|≤0.18 VIF 1.02–1.05 (seasonal offset decorrelates ENSO/IOD — the low VIF is real). (5) **negative VP direct components in ~21–27% of cells** (was silently clamped) — now reported honestly, direct read as distribution+depth-scaling not per-cell. (6) first-difference F ≥ linear-detrend F (0.33–0.37). **Still open (needs external data): observational validation (satellite LSWT / in-situ river T) and ISIMIP multi-MODEL ensemble — flagged as future work in v2 caveats.**

**FAIR+R release metadata (2026-07-14):** repo-root files `README.md` (rewritten, FAIR sections), `CITATION.cff`, `.zenodo.json`, `codemeta.json`, `LICENSE` (MIT, code), `LICENSE-CC-BY-4.0.md` (data/figs/text), `install.R` (R deps), `input/DATA_SOURCES.md` (full provenance). License choice: **MIT code + CC-BY-4.0 data** (Zenodo primary cc-by-4.0). Placeholders to fill: ORCIDs, Zenodo DOI, GitHub URL, release date. `.gitignore` now excludes large third-party input NetCDFs (clim_data/laketemp/rivertemp *.nc → ISIMIP/Zenodo, two-tier deposit) + R/TeX artifacts. Data-index citations added to §2.1: golub2022 (GMD 15,4597), frieler2024 (GMD 17,1), + NOAA data-source refs oni_cpc/nao_cpc/iod_psl (URLs from `input/clim_index/reference`). Bai/Chiew/Kraemer/Williams added earlier for intro para 2. Removed Köppen-Geiger (dead code in 04 + all zone-specific text claims: "0.74 cold lakes", "13% temperate rivers", "tropical maritime wind-mixing", climate-zone evidence).

**Annual-mean index sensitivity (2026-07-14):** `code/11_build_annual_indices.R` → `input/clim_index/annual_indices.csv` (calendar-year means, SAME column names enso_djf/nao_djfm/iod_son for engine reuse). `03`/`04` + figure `F2`/`F3` now take env-vars **FP_INDEX_FILE** + **FP_OUT_ROOT** (+ FIG_OUT for figures); defaults reproduce main analysis exactly (zero risk). Annual variant → `output/footprint_annual/` + SI Figs **S9 (footprint)/S10 (mechanism)** = `output/figures/supplementary/FS_annual_{footprint,mechanism}.png`. Result: annual indices give BROADER but WEAKER, less-ENSO-specific footprint (surface 7.4%→20.1% extent, F 0.30→0.20; ENSO dominance 78%→44%, NAO share 10%→44%) but UNCHANGED mechanism (mediated 0.99/0.93/0.98) — justifies season-specific choice. Annual-vs-seasonal index r: ENSO 0.41, NAO 0.53, IOD 0.90. Referenced in v4 §2.2 + new SI section.

**Interactive method explainer (2026-07-14):** `output/figures/interactive/footprint_method_explainer.html` — self-contained JS that recomputes F live on 4 real cells (Andes-ENSO F=0.64, Amazon 0.52, N-Atlantic-NAO 0.16, NZ-weak 0.02); reproduces pipeline exactly. Published Artifact https://claude.ai/code/artifact/8db7714a-aacf-4a1e-8e6f-afb255cebe59 . Built with Okabe-Ito ENSO/IOD/NAO colors, serif/sans/mono, light+dark. Method-outdated Q: linear-additive is fine for a variance-partitioning hypothesis (n=72, parsimony), the one real gap is asymmetry/interaction (Williams&Hanan destructive interference, Bai nonlinear) — offered but not yet built.

**Repo hygiene (2026-07-13/14):** river NetCDFs compressed (rivtemp annual 120→22MB deflate; triver monthly→annual); oversized originals moved to gitignored `_largefiles/` (NOTE: user renamed `_largefiles`→`largefiles` in script 02). All 11 scripts made path-reproducible (auto repo-root resolver replacing hardcoded `/home/dmercado/...`; walks up to find input/+code/). `input/lake_bathymetry/` was unused before revision #2. Removed stale empty `output/figures/fig3/` (06b now writes to `interactive/`).

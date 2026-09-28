# Fronts and Marine Heatwaves — GRL 2026

This folder contains the code and supporting data used to reproduce figures in the submitted paper on marine heatwaves (MHWs) and atmospheric fronts.

## Folder overview

- `code/` — analysis scripts, figure notebooks, shared startup configuration, and data-preparation documentation.
- `data/` — prepared input and intermediate data products used by the analysis.
- `plots/` — generated manuscript and supplementary figures.

## Code

| File | Purpose |
|---|---|
| `code/2026_coldfronts_mhws_startup.py` | Shared imports, region definitions, plotting defaults, coordinate/projection setup, and reusable mapping/plotting helpers. |
| `code/data_prep_manuscript.py` | Prepares frontal events and 6-hourly MHW data, extracts ERA5 front strength, and builds the MHW and frontal catalogs and frontal metrics. |
| `code/data_prep_manuscript.md` | Methods and step-by-step description of the manuscript data-preparation workflow and its outputs. |
| `code/Fig1ab_intro_map.ipynb` | Creates the introductory regional map and Figure 1 map elements. |
| `code/Fig1_and_3_stats.py` | Analyzes the MHW and frontal catalogs and creates regional/seasonal statistics and plots for Figures 1 and 3. |
| `code/Fig2_composite_stderr_figure.ipynb` | Builds event-centered frontal temperature composites and associated standard-error comparisons for Figure 2 and supplementary seasonal analyses. |
| `code/Fig4_ERA5_composites.ipynb` | Prepares ERA5 fields and plots atmospheric composites associated with frontal passages for Figure 4. |
| `code/SI_mixed_layer_budget.py` | Estimates the surface heat flux needed to produce a specified temperature change over a selected interval and mixed-layer depth; creates the supplementary heat-budget diagnostic. |

## Data

| File | Purpose |
|---|---|
| `data/ERA5_loc.nc` | ERA5 fields extracted at the study locations, including sea-level pressure, near-surface temperature, wind stress, and heat fluxes. |
| `data/M1_dsig100_rolling36h.nc` | Jordan Basin mooring stratification time series used in frontal-response calculations. |
| `data/Pioneer_BulkStratInshore_rolling36h.nc` | Pioneer mooring bulk-stratification data used in frontal-response calculations. |
| `data/MHW_days_2010_2022.nc` | Gridded MHW-day data for the period indicated in the filename, used in spatial/seasonal analysis. |
| `data/coldfront_xr.nc` | Cold-front binary time series by region and time. |
| `data/warmfront_loc.nc` | Warm-front binary time series at study locations. |
| `data/coldfront_strength_loc.nc` | Cold-front strength time series at study locations. |
| `data/frontal_freq_fig1.nc` | Frontal-frequency data used in Figure 1. |
| `data/interm_loc_6hr_oisst.npy` | Pickle-serialized intermediate 6-hourly OISST/MHW time series used by composite analyses. |
| `data/mhw_loc_6hr_oisst.npy` | Pickle-serialized MHW event records detected from 6-hourly OISST. |
| `data/front_metrics_era5_full.npy` | Pickle-serialized full-record frontal-event metrics, including event times and durations. |
| `data/front_catalog.npy` | Pickle-serialized frontal-passage catalog with strength, temperature-response, stratification, and coincident-MHW fields. |
| `data/mhw_catalog.npy` | Pickle-serialized catalog with one record per MHW event and its timing and intensity metrics. |
| `data/era5_msl_front_mean_loc.npy` | Pickle-serialized mean sea-level-pressure fields associated with cold-front passages by location. |

The `.npy` files listed above contain Python pickle-serialized objects, not plain NumPy arrays. Load them with the project’s `pickle_load` helper or `pickle.load`, rather than treating them as ordinary NumPy array files.

## Generated figures

### Main figures

- `plots/fig1b_fronts_annual.jpg` — annual frontal-passage plot.
- `plots/fig1_seasonal_fronts_heatmap.jpg` — seasonal frontal-passage heatmap.
- `plots/fig1_seasonal_mhw_days_heatmap.jpg` — seasonal MHW-day heatmap.
- `plots/fig2_composite_stderr_tempa_dt_-24_36hrs.jpg` — Figure 2 temperature composites and standard errors.
- `plots/fig3_stratification_change_seasonal_dt24_36.jpg` — seasonal stratification-change plot.
- `plots/fig3_stratification_seasonal.jpg` — seasonal stratification plot.
- `plots/fig4_era5_composites_MAB_North.jpg` — Figure 4 ERA5 composites for MAB North.

### Supplementary figures

- `plots/SI_composite_destratified_stderr_tempa_dt_-24_36hrs.jpg` — composites for the unstratified season group.
- `plots/SI_composite_stratified_stderr_tempa_dt_-24_36hrs.jpg` — composites for the stratified season group.
- `plots/SI_composite_stderr_tempa_dt_-96_96hrs.jpg` — composites using the wider lag interval.
- `plots/SI_delta_t_peak_seasonal_PI_JB_dt96_24.jpg` — seasonal peak-related temperature-change comparison for Pioneer and Jordan Basin.

## Reproducing the analysis

The notebooks and scripts expect the project’s Python environment and shared startup configuration. Run analysis code from the `code/` directory unless a notebook or script specifies otherwise, since many paths are relative.

Some original source datasets are not included here. In particular, the preparation and ERA5 workflows reference OISST and ERA5 files stored at external project/HPC paths. Re-running the entire workflow therefore requires access to those source datasets and the relevant Python packages; the included `data/` files are prepared or intermediate products used by the figure workflows.

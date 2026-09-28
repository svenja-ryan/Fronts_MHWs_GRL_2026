# Manuscript Data Preparation

## Purpose

`data_prep_manuscript.py` prepares cold-front and marine heatwave (MHW) data for the manuscript analysis. It combines frontal-event identification, 6-hourly MHW detection, frontal metrics, and event-catalog construction for eight study regions.

The script should be run from the `scripts/` directory because its input and output paths are relative to that location.

## Required Inputs

The script uses:

- `../data/coldfront_MHW.mat`: binary cold-front time series and MATLAB timestamps.
- `/vast/clidex/data/obs/SST/OISST/OISSTv2.1/oisst.day.mean_1982-2022.10.nc`: daily OISST data.
- `/vast/clidex/data/reanalysis/ERA/ERA5/subsets/atmos_eoas/F_diagnostic/front_strength/strength_*.nc`: ERA5 cold-front strength files.
- `../data/M1_dsig100_rolling36h.nc`: rolling-density stratification for Jordan Basin.
- `../data/Pioneer_BulkStratInshore_rolling36h.nc`: rolling-density stratification for Pioneer.
- `../utils/2026_coldfronts_mhws_startup.py`: shared imports and analysis configuration.

The eight locations are Jordan Basin, Scotian Shelf 1, Scotian Shelf 2, Pioneer, Slope North, Slope South, MAB North, and MAB South.

## Processing Steps

### 1. Prepare frontal events

The MATLAB frontal data are loaded and converted to an `xarray` dataset. Coordinates are assigned for the eight locations, and the record is restricted to 1983--2018. February 29 is removed so that the frontal and MHW time series have matching non-leap-day calendars.

Each contiguous cold-front occurrence is assigned an event number at each location. Event IDs greater than zero represent individual frontal passages; zero represents no front.

### 2. Interpolate SST to the frontal time grid

Daily OISST is converted from 0--360 to -180--180 longitude coordinates and sampled at the nearest grid point for each study location. The daily SST series are linearly interpolated to the frontal timestamps, producing 6-hourly SST series for 1983--2018.

### 3. Detect MHWs

MHW thresholds and events are calculated with `xmhw` using:

- Climatology period: 1983--2012
- Time step: 6-hourly
- Smoothing width: 121 time steps, equivalent to approximately 30 days
- Maximum gap within an event: 8 time steps
- Minimum duration: 20 time steps, equivalent to 5 days

The intermediate MHW output is augmented with temperature anomaly relative to the seasonal climatology. These anomalies are used later for temperature-change calculations around frontal passages.

### 4. Calculate frontal metrics

For every numbered front at each location, the script records:

- Event ID
- Start time
- End time
- Duration in hours

The full-record metric dictionary is written as `front_metrics_era5_full.npy` when the corresponding `pickle_save` line is enabled in the script.

### 5. Extract frontal strength

The ERA5 frontal-strength files are opened and the nearest grid point is selected for each study location. The resulting data are combined into `front_strength_loc`, restricted to 1983--2018, and filtered to remove February 29. The per-location data are stored in `front_strength` and used directly when constructing the frontal catalog. The script does not require the previously generated `coldfront_strength_loc.nc` file.

### 6. Build the MHW catalog

One row is created for each detected MHW event. The catalog includes:

- Region and event ID
- Start, end, and peak dates
- Duration in days
- Maximum and mean intensity
- Cumulative intensity in degree-days
- Meteorological season based on the peak month

### 7. Build the frontal catalog

One row is created for each valid frontal passage. The catalog includes frontal timing, duration, season, front strength at onset, maximum strength during the event, and the time of maximum strength.

For each passage, the script also calculates:

- Temperature anomaly change from 24 hours before to 36 hours after the front
- Temperature anomaly change from 5 days before to 5 days after the front
- Temperature anomaly change from 5 days before to 1 day before the front
- SST before and after the front
- Stratification before, during, and after the front where mooring data are available

MHW metadata are attached when an MHW is present during the 24 hours before the frontal passage. MHW event type is classified as:

- `short`: duration less than 10 days
- `long`: duration of at least 10 days

Finally, frontal strength is categorized separately within each region using percentile thresholds:

- `weak`: 0th--25th percentile
- `moderate`: 25th--75th percentile
- `strong`: 75th--90th percentile
- `extreme`: 90th--100th percentile

## Outputs

The script writes the following files to `../data/` using the repository's `pickle_save` helper. Although the files use the `.npy` extension, they contain Python pickle-serialized objects:

- `mhw_catalog.npy`: pandas DataFrame with one row per MHW event.
- `front_catalog.npy`: pandas DataFrame with one row per valid frontal passage, with coincident MHW fields where applicable.
- `front_metrics_era5_full.npy`: dictionary of frontal-event metrics keyed by region.

Because these files are Python-specific serialized objects, they should be loaded with the repository's `pickle_load` helper or `pickle.load` rather than `numpy.load` alone.

Example:

```python
from utils.general import pickle_load

mhw_catalog = pickle_load("mhw_catalog", "../data/")
front_catalog = pickle_load("front_catalog", "../data/")
front_metrics = pickle_load("front_metrics_era5_full", "../data/")
```

# %% 
# MHW & Cold Front Statistics
#
# This script prepares the frontal-event and 6-hourly MHW data, computes
# full-record frontal metrics, and builds manuscript-ready event catalogs.
#
# Outputs are saved as pickled files:
#   ../data/mhw_catalog.npy
#   ../data/front_catalog.npy
#   ../data/front_metrics_era5_full.npy

import os, sys
import scipy.io as sio

STARTUP_2026_coldfronts_mhws = "../utils/2026_coldfronts_mhws_startup.py"
exec(open(STARTUP_2026_coldfronts_mhws).read())

sys.path.append("../")
from utils.datafun import matlab2datetime
from utils.general import pickle_save
from xmhw.xmhw import threshold, detect


# =============================================================================
# LOAD AND PREPARE FRONTAL DATA
# =============================================================================
## read frontal data from Rhys' file
_ds = sio.loadmat('../data/coldfront_MHW.mat', squeeze_me=True)

# closest ERA 5 locations to initial locations --> provided by Rhys
lat_fronts = [43.5,43,44.5,40.25,41.5,39,40.75,38.5]
lon_fronts = np.array([292,295,296.5,289,295,289,288,285.5])-360
fronts = _ds['coldfronts']
time = _ds['timenum']

# convert time to datetime
timefront=[]
for i in range(len(time)):
    timefront.append(matlab2datetime(time[i]))

# create xarray
front = xr.DataArray(data=fronts,
                  dims=["loc","time"],
                  coords=dict(lon=(["loc"], lon_fronts),
                              lat=(["loc"], lat_fronts),
                              loc=(["loc"],['JordanBasin','ScotianShelf1','ScotianShelf2','Pioneer','SlopeN','SlopeS','MAB_North','MAB_south']),
                              time=timefront)
                 )

# remove Feb 29 to make compatible with MHW data
front = front.sel(time=~((front.time.dt.month == 2) &
                                (front.time.dt.day == 29))).sel(time=slice('1983','2018'))

# delete variables not used
del timefront,time,_ds,lat_fronts,lon_fronts

## identifiy separate events
front = front.to_dataset(name='front')
frontnr = front['front'].astype(int).values
for i in range(8): # loop through locations
    # select one region
    dummy = front['front'][i,:].astype(int).values
    # find start and end of each event
    wo_start =  np.where(np.diff(dummy)==1)[0]+1
    wo_end =  np.where(np.diff(dummy)==-1)[0]+1 # have to add +1 due non-inclusive selection

    # add new vector with event numbering
    for j in range(len(wo_start)):
        frontnr[i,wo_start[j]:wo_end[j]]=j+1

# add to dataset
front['event'] = (('loc','time'),frontnr)

# delete unused variables
del frontnr,wo_start,wo_end




#%% =============================================================================
# LOAD AND PREPARE MHW DATA
# =============================================================================

# extract locations in MHW dataset
names = ['JordanBasin','ScotianShelf1','ScotianShelf2','Pioneer','SlopeN','SlopeS','MAB_North','MAB_south']
lat = [43.49,43,44.5,40.362,41.5,39,40.693,38.5]
lon = [-67.88,-65,-63.4,-70.8785,-65,-70.8785,-72.049,-74.5]


# interpolate OISST to 6 hourly and apply MHW detection
oisst = xr.open_dataset('/vast/clidex/data/obs/SST/OISST/OISSTv2.1/oisst.day.mean_1982-2022.10.nc')
oisst['lon'] = oisst['lon']-360

oisst_loc = {}
loc={}
for name,i in zip(names,range(len(names))):
    loc[name] = [lat[i],lon[i]]
    oisst_loc[name] = oisst.sel(lon=lon[i],lat=lat[i],method='nearest').load()

# interpolate from daily to 6 hourly to match frontal data for time period of frontal data (1983 - 2018)
oisst_loc_6hr = {}
for name,i in zip(names,range(len(names))):
    oisst_loc_6hr[name] = oisst_loc[name]['sst'].interp(time=front.time).sel(time=~((front.time.dt.month == 2) &
                                (front.time.dt.day == 29))).sel(time=slice('1983','2018'))

# MHW detection
clim_loc_6hr = {}
mhw_loc_6hr = {}
interm_loc_6hr = {}
for name,i in zip(names,range(len(names))):
    # adapted smoothing window so that it is 30days with 6 hourly timestep
    clim_loc_6hr[name] = threshold(oisst_loc_6hr[name],climatologyPeriod=['1983','2012'],
                                   smoothPercentileWidth=(30*4)+1, tstep=True)
    mhw_loc_6hr[name], interm_loc_6hr[name] = detect(oisst_loc_6hr[name].drop(('lat','lon')),
                                                     clim_loc_6hr[name]['thresh'],
                                                     clim_loc_6hr[name]['seas'],
                                                     intermediate=True,maxGap=8, tstep=True,minDuration=5*4)
    interm_loc_6hr[name] = interm_loc_6hr[name].rename({"index": "time"})

# add temperature anomaly by subtracting daily climatology
for name in names:
    dummy = interm_loc_6hr[name].copy(deep=True)
    ts_times = dummy['ts'].time
    # tile across full record
    doy = np.tile(clim_loc_6hr[name]['doy'].values, len(ts_times) // 1460) #//makes integer
    dummy = dummy.assign_coords(doy=("time", doy))
    interm_loc_6hr[name]['tempa'] = (dummy['ts'].groupby("doy") - clim_loc_6hr[name]['seas']).drop(('lat','lon')).load()

# save data needed by the catalogs
# pickle_save('interm_loc_6hr_oisst','../data/',interm_loc_6hr)
# pickle_save('clim_loc_6hr_oisst','../data/',clim_loc_6hr)
# pickle_save('mhw_loc_6hr_oisst','../data/',mhw_loc_6hr)









# %% =============================================================================
#CREATE FULL-RECORD FRONTAL METRICS to be used for composites in Figure 2
# =============================================================================
# ---------------------------------------------------------
# Create frontal metrics dictionary (similar structure to MHW metrics)
# ---------------------------------------------------------
def front_metrics(ds):
    """
    Compute basic metrics for frontal events at each location.
    """
    front_metrics_loc = {}

    for name in names:
        # Extract unique event IDs at this location
        # (ignore zeros = no event)
        dummy = np.unique(ds['event'].sel(loc=name))

        metric = {}
        metric['event'] = dummy[dummy > 0]

        # allocate lists for metrics
        start = []
        end = []
        duration = []

        # Loop over each event ID
        for i in range(len(metric['event'])):
            dummy = front.sel(loc=name) # select data for this location
            event_data = dummy.where(dummy['event'] == metric['event'][i], drop=True) # isolate current event
            start.append(event_data['time'].min().values) # start time = first occurrence
            end.append(event_data['time'].max().values) # end time = last occurrence
            duration.append(len(event_data['time']) * 6)  # duration = number of time steps * 6 hours (assuming 6-hourly data)

        # store metrics
        metric['time_start'] = start
        metric['time_end'] = end
        metric['duration'] = duration
        front_metrics_loc[name] = metric

    return front_metrics_loc

front_metrics_full = front_metrics(front)
# pickle_save('front_metrics_era5_full','../data/',front_metrics_full)



#%% =============================================================================
# Extract frontal strength data DATA USED BY THE FRONT CATALOG
# =============================================================================

front_strength_data = xr.open_mfdataset(
    "/vast/clidex/data/reanalysis/ERA/ERA5/subsets/atmos_eoas/F_diagnostic/front_strength/strength_*.nc",combine="by_coords")


lats = [43.5,43,44.5,40.25,41.5,39,40.75,38.5]
lons = np.array([292,295,296.5,289,295,289,288,285.5])-360

series = []
for name, lon0, lat0 in zip(names, lons, lats):
    da_i = front_strength_data["F_strength"].sel(
        time=slice("1983-01-01", "2018-12-31")).sel(
        lon=lon0 % 360, # change to 0-360 for ERA5 data
        lat=lat0,
        method="nearest",
    ).compute()

    da_i = da_i.expand_dims(loc=[name]) # add location dimension

    # assign coordinates for lon and lat
    da_i = da_i.assign_coords(
        lon=("loc", [lon0]),
        lat=("loc", [lat0]),
    )

    # append extracted data to list
    series.append(da_i)

front_strength_loc = xr.concat(series, dim="loc")
front_strength_loc = front_strength_loc.transpose("loc", "time")

front_strength_loc.attrs["comment"] = "Extracted from Rhys front strength data at specified locations. Data at /vast/clidex/data/reanalysis/ERA/ERA5/subsets/atmos_eoas/F_diagnostic/front_strength/"

# remove Feb 29 to make compatible with MHW data
front_strength_loc = front_strength_loc.sel(time=~((front_strength_loc.time.dt.month == 2) &
                                (front_strength_loc.time.dt.day == 29)))

# front_strength_loc.to_netcdf('../data/coldfront_strength_loc.nc')
front_strength = {name: front_strength_loc.sel(loc=name) for name in names}






# %%=============================================================================
#  LOAD Mooring data DATA USED BY THE FRONT CATALOG
# =============================================================================
rho_JB = xr.open_dataset('../data/M1_dsig100_rolling36h.nc')['sigma_t']
rho_PI = xr.open_dataset('../data/Pioneer_BulkStratInshore_rolling36h.nc')['ddens']


def _season_label(month):
    if month in [12, 1, 2]:  return 'DJF'
    if month in [3, 4, 5]:  return 'MAM'
    if month in [6, 7, 8]:  return 'JJA'
    return 'SON'



#%% =============================================================================
# BUILD MHW CATALOG
# =============================================================================
records = []
for name in names:
    ds = mhw_loc_6hr[name]
    interm = interm_loc_6hr[name]
    n = len(ds['event'])

    ts_start_vals = ds['time_start'].values
    ts_end_vals = ds['time_end'].values
    ts_peak_vals = ds['time_peak'].values

    for i in range(n):
        t_start = pd.Timestamp(ts_start_vals[i])
        t_end = pd.Timestamp(ts_end_vals[i])
        t_peak = pd.Timestamp(ts_peak_vals[i])

        records.append({
            'region':               name,
            'event':                int(ds['event'].values[i]),
            'date_start':           t_start,
            'date_end':             t_end,
            'date_peak':            t_peak,
            'duration':             float(ds['duration'].values[i]) / 4,
            'intensity_max':        float(ds['intensity_max'].values[i]),
            'intensity_mean':       float(ds['intensity_mean'].values[i]),
            'intensity_cumulative': float(ds['intensity_cumulative'].values[i]) / 4,
            'season':               _season_label(t_peak.month),
        })

mhw_catalog = pd.DataFrame(records)
# pickle_save('mhw_catalog', '../data/', mhw_catalog)






#%% =============================================================================
# BUILD FRONT CATALOG
# =============================================================================
records_front = []

for name in names:
    ds_mhw = mhw_loc_6hr[name]
    interm = interm_loc_6hr[name]
    dummy_front = front_metrics_full[name]
    dummy_strength = front_strength[name]

    front_starts_all = pd.to_datetime(dummy_front['time_start'])
    strength_by_pos = dummy_strength.sel(
        time=xr.DataArray(np.asarray(front_starts_all), dims='e'),
        method='nearest',
    ).values.astype(float)

    n_event_tmp = len(dummy_front['event'])
    strength_max_by_pos = np.full(n_event_tmp, np.nan)
    strength_peak_time_by_pos = np.full(n_event_tmp, np.datetime64('NaT'), dtype='datetime64[ns]')
    for i in range(n_event_tmp):
        t0 = dummy_front['time_start'][i]
        t1 = dummy_front['time_end'][i]
        strength_window = dummy_strength.sel(time=slice(t0, t1))
        if strength_window.size:
            strength_max_by_pos[i] = float(strength_window.max())
            peak_idx = int(strength_window.argmax())
            strength_peak_time_by_pos[i] = strength_window['time'].values[peak_idx]

    event_ids_all = ds_mhw['event'].values.astype(int)
    mhw_meta = {}
    for idx, eid in enumerate(event_ids_all):
        t_start = pd.Timestamp(ds_mhw['time_start'].values[idx])
        t_end = pd.Timestamp(ds_mhw['time_end'].values[idx])
        t_peak = pd.Timestamp(ds_mhw['time_peak'].values[idx])
        mhw_meta[eid] = {
            'mhw_event': eid,
            'mhw_date_start': t_start,
            'mhw_date_end': t_end,
            'mhw_date_peak': t_peak,
            'duration': float(ds_mhw['duration'].values[idx]) / 4,
            'intensity_max': float(ds_mhw['intensity_max'].values[idx]),
            'intensity_mean': float(ds_mhw['intensity_mean'].values[idx]),
            'intensity_cumulative': float(ds_mhw['intensity_cumulative'].values[idx]) / 4,
            'mhw_season': _season_label(t_peak.month),
            'mhw_type': 'short' if float(ds_mhw['duration'].values[idx]) / 4 < 10 else 'long',
        }

    n_time = 41
    n_event = len(dummy_front['event'])
    mhw_event_comp = np.full((n_time, n_event), np.nan)
    valid_event = np.full(n_event, False)

    for i in range(n_event):
        t0 = dummy_front['time_start'][i] - np.timedelta64(5, 'D')
        t1 = dummy_front['time_start'][i] + np.timedelta64(5, 'D')
        subset = interm.sel(time=slice(t0, t1))
        if len(subset['tempa']) == n_time:
            mhw_event_comp[:, i] = subset['events'].values
            valid_event[i] = True

    for i in range(n_event):
        if not valid_event[i]:
            continue

        front_t_start = pd.Timestamp(dummy_front['time_start'][i])
        front_t_end = pd.Timestamp(dummy_front['time_end'][i])
        front_dur_h = dummy_front['duration'][i]
        front_season = _season_label(front_t_start.month)

        t_pre = front_t_start - pd.Timedelta(hours=24)
        t_post = front_t_start + pd.Timedelta(hours=36)

        rho_pre = rho_post = rho_front = delta_rho = np.nan
        try:
            ta_pre = float(interm['tempa'].sel(time=t_pre, method='nearest'))
            ta_post = float(interm['tempa'].sel(time=t_post, method='nearest'))
            sst_pre = float(interm['ts'].sel(time=t_pre, method='nearest'))
            sst_post = float(interm['ts'].sel(time=t_post, method='nearest'))
            delta_t_24_36 = ta_post - ta_pre
            delta_t_96_96 = (float(interm['tempa'].sel(time=front_t_start + pd.Timedelta(days=5), method='nearest')) -
                              float(interm['tempa'].sel(time=front_t_start - pd.Timedelta(days=5), method='nearest')))
            delta_t_peak_96_24 = (float(interm['tempa'].sel(time=front_t_start - pd.Timedelta(days=1), method='nearest')) -
                                   float(interm['tempa'].sel(time=front_t_start - pd.Timedelta(days=5), method='nearest')))
            if name == "JordanBasin":
                rho_pre = float(rho_JB.sel(time=t_pre, method='nearest'))
                rho_post = float(rho_JB.sel(time=t_post, method='nearest'))
                delta_rho = rho_post - rho_pre
                rho_front = float(rho_JB.sel(time=front_t_start, method='nearest'))
            elif name == "Pioneer":
                rho_pre = float(rho_PI.sel(time=t_pre, method='nearest'))
                rho_post = float(rho_PI.sel(time=t_post, method='nearest'))
                delta_rho = rho_post - rho_pre
                rho_front = float(rho_PI.sel(time=front_t_start, method='nearest'))
        except Exception:
            ta_pre = ta_post = sst_pre = sst_post = delta_t_24_36 = delta_t_96_96 = delta_rho = rho_pre = rho_post = delta_t_peak_96_24 = rho_front = np.nan

        eids = np.unique(mhw_event_comp[16:21, i])
        eids = eids[~np.isnan(eids)].astype(int)

        base_row = {
            'region': name,
            'front_event': int(dummy_front['event'][i]),
            'front_start': front_t_start,
            'front_end': front_t_end,
            'front_duration': front_dur_h,
            'front_season': front_season,
            'front_strength': strength_by_pos[i],
            'front_strength_max': strength_max_by_pos[i],
            'front_peak_time': strength_peak_time_by_pos[i],
            'delta_t_24_36': delta_t_24_36,
            'delta_t_96_96': delta_t_96_96,
            'delta_t_peak_96_24': delta_t_peak_96_24,
            'sst_pre': sst_pre,
            'sst_post': sst_post,
            'delta_rho': delta_rho,
            'rho_pre': rho_pre,
            'rho_post': rho_post,
            'rho_front': rho_front,
        }

        if len(eids) == 0:
            records_front.append(base_row)
        else:
            for eid in eids:
                if eid not in mhw_meta:
                    continue
                row = base_row.copy()
                row.update(mhw_meta[eid])
                records_front.append(row)

front_catalog = pd.DataFrame(records_front)

# Define labels for frontal strength
labels = ["weak", "moderate", "strong", "extreme"]

# Assign percentile-based categories within each region.
front_catalog["front_type"] = (
    front_catalog
    .groupby("region")["front_strength"]
    .transform(lambda x: pd.qcut(x, q=[0, 0.25, 0.75, 0.90, 1.0],
                                 labels=labels))
)

# Move the new category to directly follow the strength column.
col = front_catalog.pop("front_type")
idx = front_catalog.columns.get_loc("front_strength") + 1
front_catalog.insert(idx, "front_type", col)

# pickle_save('front_catalog', '../data/', front_catalog)

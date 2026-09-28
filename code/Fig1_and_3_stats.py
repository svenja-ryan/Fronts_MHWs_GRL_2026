#%% =============================================================================
# FIGURES 1 AND 3: MHW AND COLD-FRONT STATISTICS
# =============================================================================
#
# Overview
# --------
# This script loads the MHW and frontal-passage catalogs, computes regional
# summary statistics, and creates manuscript figures for frontal temperature
# response, stratification changes, seasonal MHW days, and frontal passages.
#
# The startup file provides the shared scientific Python imports, plotting
# helpers, region names, and figure styling used throughout this script.
#
# Main outputs
# ------------
# - Delta-T and stratification whisker plots
# - Seasonal comparison plots for MHW states
# - Monthly heatmaps of MHW days and frontal passages
#
#
# =============================================================================
# SETUP AND DATA LOADING
# =============================================================================
import os, sys
STARTUP_2026_coldfronts_mhws = "./2026_coldfronts_mhws_startup.py"
exec(open(STARTUP_2026_coldfronts_mhws).read())
from matplotlib.lines import Line2D
import seaborn as sns
import cmocean as cmo


# Load prepared catalogs
front_catalog = pickle_load('front_catalog', '../data/')
mhw_catalog = pickle_load('mhw_catalog', '../data/')


# =============================================================================
# DELTA-T STATISTICS FOR FRONTAL PASSAGES
# =============================================================================

# -----------------------------------------------------------------------------
# Summary statistics: delta_t by region and MHW type
# -----------------------------------------------------------------------------
def var_stats(ds=front_catalog,var=None):
    rows = []
    for name in ds['region'].unique():
        sub = ds[ds['region'] == name]
        groups = [
            ('No MHW', sub[sub['mhw_event'].isna()][var].dropna()),
            ('Short',  sub[sub['mhw_type'] == 'short'][var].dropna()),
            ('Long',   sub[sub['mhw_type'] == 'long'][var].dropna()),
        ]
        for label, grp in groups:
            rows.append({'Region': name, 'MHW type': label,
                        'n': len(grp), 'mean': grp.mean(), 'median': grp.median(), 'std': grp.std()})

    return pd.DataFrame(rows).set_index(['Region', 'MHW type']).round(3)


# -----------------------------------------------------------------------------
# Plot all MHW states for each region
# -----------------------------------------------------------------------------
def plot_var_stats_box_whisker_all(var_stats_df,var_name,names=names):
    """
    Box-and-whisker plot of a variable (delta_t or rho_pre) by region and MHW type.

    Parameters
    ----------
    var_stats_df : DataFrame with index ['Region', 'MHW type'] and columns ['n', 'mean', 'median', 'std']
    var_name : str — name of the variable for labeling
    names : list of region names
    plotsave : str — filename to save the plot (optional)
    """
    mhw_states = ["No MHW", "Short", "Long"]
    colors = {"No MHW": "#4C72B0", "Short": "#DD8452", "Long": "#C44E52"}
    offsets = {"No MHW": 0.22, "Short": 0.0, "Long": -0.22}

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.axvline(0, color="0.3", lw=0.8, ls="--", zorder=0)

    for i, region in enumerate(names):
        for state in mhw_states:
            row = var_stats_df.loc[(region, state)]
            mean, std, n = row["mean"], row["std"], row["n"]
            se = std / np.sqrt(n)
            ax.errorbar(
                mean, i + offsets[state], xerr=se,
                fmt="o", color=colors[state],
                markersize=7, markeredgecolor="white", markeredgewidth=0.8,
                elinewidth=1.6, capsize=3, capthick=1.2, zorder=3,
            )
            # TO ANNOTATE COUNTS FOR EACH DISTRIBUTION UNCOMMENT
            if state == 'Short':  
                ax.text(mean -0.01, i + offsets[state] + 0.13, f"n={int(n)}",
                va="center", fontsize=7, color=colors[state], alpha=0.8)
            else:
                ax.text(mean -0.01, i + offsets[state] - 0.13, f"n={int(n)}",
                va="center", fontsize=7, color=colors[state], alpha=0.8)

    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    ax.set_xlabel(rf"${var_name}$ across frontal passage (°C)")
    ax.set_title("SST Anomaly response to frontal passage by MHW state")
    ax.grid(axis="x", alpha=0.25)
    ax.set_axisbelow(True)

    legend_elems = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=colors[s],
            markersize=8, markeredgecolor="white", label=s)
        for s in mhw_states
    ]
    ax.legend(handles=legend_elems, loc="upper left", frameon=True)

    fig.tight_layout()
    # plt.show()

    return fig,ax

# -----------------------------------------------------------------------------
# Plot a compact version of the MHW-state whisker figure
# -----------------------------------------------------------------------------
def plot_var_stats_box_whisker2(var_stats_df,var_name,names=names):
    """
    Box-and-whisker plot of a variable (delta_t or rho_pre) by region and MHW type.

    Parameters
    ----------
    var_stats_df : DataFrame with index ['Region', 'MHW type'] and columns ['n', 'mean', 'median', 'std']
    var_name : str — name of the variable for labeling
    names : list of region names
    plotsave : str — filename to save the plot (optional)
    """
    mhw_states = ["No MHW", "Short", "Long"]
    colors = {"No MHW": "#4C72B0", "Short": "#DD8452", "Long": "#C44E52"}
    offsets = {"No MHW": 0.22, "Short": 0.0, "Long": -0.22}

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.axhline(0, color="0.3", lw=0.8, ls="--", zorder=0)

    for i, region in enumerate(names):
        for state in mhw_states:
            row = var_stats_df.loc[(region, state)]
            mean, std, n = row["mean"], row["std"], row["n"]
            se = std /np.sqrt(n)

            ax.errorbar(
                i + offsets[state], mean, yerr=se,
                fmt="o",
                color=colors[state],
                markersize=7,
                markeredgecolor="white",
                markeredgewidth=0.8,
                elinewidth=1.6,
                capsize=3,
                capthick=1.2,
                zorder=3,
            )

    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha="right")
    # ax.set_ylabel(rf"bulk stratification at frontal passage (kg/m³)")
    # ax.set_title("Bulk stratification at frontal passage by MHW state")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    # ax.set_ylim(0,4)

    legend_elems = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=colors[s],
            markersize=8, markeredgecolor="white", label=s)
        for s in mhw_states
    ]
    ax.legend(handles=legend_elems, loc="lower right", frameon=True)

    fig.tight_layout()
    # plt.show()
    return fig,ax

# -----------------------------------------------------------------------------
# Compare MHW states across two seasonal groups
# -----------------------------------------------------------------------------
def plot_seasonal_var_stats(var, regions, var_name):
    """Plot mean +/- SE for each MHW state and region in two season groups."""
    season_groups = {
        "DJF + MAM": ["DJF", "MAM"],
        "JJA + SON": ["JJA", "SON"],
    }
    mhw_states = ["No MHW", "Short", "Long"]
    state_colors = {"No MHW": "#4C72B0", "Short": "#DD8452", "Long": "#C44E52"}
    season_markers = {"DJF + MAM": "o", "JJA + SON": "s"}
    season_alpha = {"DJF + MAM": 0.5, "JJA + SON": 1}
    state_offsets = {"No MHW": -0.24, "Short": 0.0, "Long": 0.24}
    season_offsets = {"DJF + MAM": -0.07, "JJA + SON": 0.07}

    fig, ax = plt.subplots(figsize=(5, 3))
    ax.axhline(0, color="0.3", lw=0.8, ls="--", zorder=0)

    for region_index, region in enumerate(regions):
        for state in mhw_states:
            for season_label, season_names in season_groups.items():
                seasonal_df = front_catalog[
                    front_catalog["front_season"].isin(season_names)
                ]
                stats = var_stats(ds=seasonal_df, var=var)

                if (region, state) not in stats.index:
                    continue

                row = stats.loc[(region, state)]
                n = row["n"]
                if n == 0:
                    continue
                se = row["std"] / np.sqrt(n)
                x_position = (
                    region_index
                    + state_offsets[state]
                    + season_offsets[season_label]
                )
                ax.errorbar(x_position, row["mean"], yerr=se,
                            fmt=season_markers[season_label], color=state_colors[state],
                            markersize=6, markeredgecolor="white", markeredgewidth=0.8,
                            elinewidth=1.4, capsize=3, capthick=1.1, zorder=3,
                            alpha=season_alpha[season_label])

    ax.set_xticks(range(len(regions)))
    ax.set_xticklabels(regions)#, rotation=45, ha="right")
    ax.set_ylabel(rf"${var_name}$ ")
    # ax.set_title("Seasonal response to frontal passage by MHW state")
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)

    state_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=state_colors[state],
               markersize=7, markeredgecolor="white", label=state)
        for state in mhw_states
    ]
    season_handles = [
        Line2D([0], [0], marker=season_markers[season_label], color="0.3", linestyle="None",
               markersize=6, label=season_label)
        for season_label in season_groups
    ]
    ax.legend(handles=state_handles + season_handles, loc="best", frameon=True, ncol=2,fontsize=6)
    fig.tight_layout()
    return fig, ax






# =============================================================================
#%% CREATE DELTA-T FIGURES
# =============================================================================
# font_for_pres()
#-------------------------------------------------------------------------------
# create stats for Delta T
# can subset catalog to only fronts of specific types if desired, e.g., weak fronts
# df_sub = front_catalog#[front_catalog["front_type"].isin(["strong","extreme"])].copy()
df_sub = front_catalog#[front_catalog["front_season"].isin(season)].copy()
delta_t_stats = var_stats(ds=df_sub,var='delta_t_24_36')
delta_t_stats2= var_stats(ds=df_sub,var='delta_t_96_96')
delta_t_stats3 = var_stats(ds=df_sub,var='delta_t_peak_96_24')
    
# plot box-and-whisker for delta_t
# fig,ax = plot_var_stats_box_whisker2(delta_t_stats,var_name=r"\Delta T_{24-36}",names=['Pioneer','JordanBasin','MAB_North'])
# fig,ax = plot_var_stats_box_whisker2(delta_t_stats2,var_name=r"\Delta T_{96-96}",names=['Pioneer','JordanBasin','MAB_North'])
fig,ax = plot_var_stats_box_whisker2(delta_t_stats3,var_name=r"\Delta T_{peak_96-24}",names=['Pioneer','JordanBasin','MAB_North'])
ax.set_ylabel(rf"$\Delta_\tau T$ [°C]")


#-------------------------------------------------------------------------------
# Seasonal comparison: DJF + MAM versus JJA + SON.
fig_seasonal, ax_seasonal = plot_seasonal_var_stats(
    var="delta_t_peak_96_24",
    regions=["Pioneer", "JordanBasin"],
    var_name=r"\Delta_\tau T \, [°C]",
)
finished_plot(fig_seasonal, "../plots/manuscript/SI_delta_t_peak_seasonal_PI_JB_dt96_24.jpg")






# =============================================================================
#%% DELTA-RHO STATISTICS FOR FRONTAL PASSAGES
# =============================================================================
# can subset catalog to only fronts of specific types if desired, e.g., weak fronts
# df_sub = front_catalog#[front_catalog["front_type"].isin(["weak"])].copy()
# rho_pre_stats = var_stats(ds=df_sub,var='rho_pre').loc[['Pioneer','JordanBasin']]
# fig,ax = plot_var_stats_box_whisker2(rho_pre_stats,var_name=rf"\Delta \rho_{{24-36}}",names=['Pioneer','JordanBasin'])
# ax.set_ylabel(rf"Stratification at frontal passage (kg/m³)")

# finished_plot(fig,ax,plotsave='../figures/fig3_bulk_stratification.jpg',dpi=300)

fig_seasonal, ax_seasonal = plot_seasonal_var_stats(
    var="rho_front",
    regions=["Pioneer", "JordanBasin"],
    var_name=r"\Delta_z\, \rho \, [kg/m³]",
)

finished_plot(fig_seasonal, "../plots/manuscript/fig3_stratification_seasonal.jpg")



# =============================================================================
# %% STRATIFICATION CHANGE DURING PASSAGE
# =============================================================================
# can subset catalog to only fronts of specific types if desired, e.g., weak fronts
# df_sub = front_catalog#[front_catalog["front_season"].isin(['JJA','SON'])].copy()
# delta_rho_stats = var_stats(ds=df_sub,var='delta_rho').loc[['Pioneer','JordanBasin']]
# fig,ax = plot_var_stats_box_whisker2(delta_rho_stats,var_name=rf"\Delta \rho_{{24-36}}",names=['Pioneer','JordanBasin'])
# ax.set_ylabel(rf"Stratification change at frontal passage (kg/m³)")
# ax.set_title(f'Season: JJA & SON')

# # 
# df_sub = front_catalog[front_catalog["front_season"].isin(['DJF','MAM'])].copy()
# delta_rho_stats = var_stats(ds=df_sub,var='delta_rho').loc[['Pioneer','JordanBasin']]
# fig,ax = plot_var_stats_box_whisker2(delta_rho_stats,var_name=rf"\Delta \rho_{{24-36}}",names=['Pioneer','JordanBasin'])
# ax.set_ylabel(rf"Stratification change at frontal passage (kg/m³)")
# ax.set_title(f'Season: DJF & MAM')
# finished_plot(fig,ax,plotsave='../figures/fig3_bulk_stratification.jpg',dpi=300)


# seasonal comparison for stratification change
fig_seasonal, ax_seasonal = plot_seasonal_var_stats(
    var="delta_rho",
    regions=["Pioneer", "JordanBasin"],
    var_name=r"\Delta_\tau\Delta_z\, \rho \, [kg/m³]",
)
finished_plot(fig_seasonal, "../plots/manuscript/fig3_stratification_change_seasonal_dt24_36.jpg")





# =============================================================================
#%% Figure 1 SEASONAL HEATMAPS: MHW DAYS AND FRONTAL PASSAGES
# =============================================================================

# ----------------------------------------------------------------------------
# Total MHW days per month for each region, accounting for event duration.
# ----------------------------------------------------------------------------
font_for_print()
region_col = "region"
selected_regions = ["JordanBasin", "Pioneer", "MAB_North"]
region_display_names = {
    "JordanBasin": "JB",
    "Pioneer": "PI",
    "MAB_North": "MAB",
}


mhw_tbl = mhw_catalog[[region_col, "date_start"]].copy()
mhw_tbl["duration"] = pd.to_numeric(mhw_catalog["duration"], errors="coerce")
mhw_tbl["date_end"] = pd.to_datetime(mhw_catalog["date_end"], errors="coerce")
mhw_tbl["date_start"] = pd.to_datetime(mhw_tbl["date_start"], errors="coerce")
mhw_tbl = mhw_tbl.dropna(subset=[region_col, "date_start"])
mhw_tbl = mhw_tbl[mhw_tbl[region_col].isin(selected_regions)].copy()

month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

alloc_rows = []
for _, row in mhw_tbl.iterrows():
    t_start = row["date_start"]
    duration_days = row.get("duration", np.nan)
    t_end = row.get("date_end", pd.NaT)

    if pd.notna(duration_days) and duration_days > 0:
        t_stop = t_start + pd.to_timedelta(float(duration_days), unit="D")
    elif pd.notna(t_end):
        # date_end is usually the final 6-hour timestep; add 6h so endpoint contributes.
        t_stop = t_end + pd.Timedelta(hours=6)
    else:
        continue

    if t_stop <= t_start:
        continue

    t_cursor = t_start
    while t_cursor < t_stop:
        month_start = pd.Timestamp(year=t_cursor.year, month=t_cursor.month, day=1)
        next_month_start = month_start + pd.offsets.MonthBegin(1)
        seg_stop = min(t_stop, next_month_start)
        seg_days = (seg_stop - t_cursor).total_seconds() / 86400.0

        alloc_rows.append(
            {
                "region": row[region_col],
                "year": t_cursor.year,
                "month": t_cursor.month,
                "mhw_days": seg_days,
            }
        )
        t_cursor = seg_stop

alloc_df = pd.DataFrame(alloc_rows)
if alloc_df.empty:
    raise ValueError("No valid MHW durations could be allocated to calendar months.")

days_by_region_month = (
    alloc_df.groupby(["region", "month"])["mhw_days"].sum().unstack(fill_value=0)
)
days_by_region_month = days_by_region_month.reindex(columns=range(1, 13), fill_value=0)
days_by_region_month = days_by_region_month.reindex(selected_regions).dropna(how="all")

if days_by_region_month.empty:
    raise ValueError("No MHW-day allocations found for selected_regions.")

heat_df = days_by_region_month.copy()
heat_df.index = heat_df.index.map(region_display_names)
heat_df.columns = month_names

start_year = int(mhw_tbl["date_start"].dt.year.min())
end_year = int(mhw_tbl["date_start"].dt.year.max())

fig, ax = plt.subplots(figsize=(4,1.5))
sns.heatmap(heat_df, ax=ax,cmap="gray_r",annot=True,fmt=".0f",linewidths=0.5,
    linecolor="white",cbar=False,
)
ax.set_title(f"Total MHW days by month ({start_year}-{end_year})")
ax.set_xlabel("")
ax.set_ylabel("")
fig.tight_layout()

finished_plot(fig,"../plots/manuscript/fig1_seasonal_mhw_days_heatmap.jpg")






# ----------------------------------------------------------------------------
# Monthly frontal-passage counts by region
# ----------------------------------------------------------------------------

font_for_print()
front_tbl = front_catalog[[region_col, "front_start"]].copy()
front_tbl["front_start"] = pd.to_datetime(front_tbl["front_start"], errors="coerce")
front_tbl = front_tbl.dropna(subset=[region_col, "front_start"])
front_tbl = front_tbl[front_tbl[region_col].isin(selected_regions)].copy()


front_by_region_month = (
    front_tbl.groupby([region_col, front_tbl["front_start"].dt.month])
    .size()
    .unstack(fill_value=0)
)
front_by_region_month = front_by_region_month.reindex(columns=range(1, 13), fill_value=0)
front_by_region_month = front_by_region_month.reindex(selected_regions).dropna(how="all")

if front_by_region_month.empty:
    raise ValueError("No monthly frontal-passage counts found for selected_regions.")

heat_front = front_by_region_month.copy()
heat_front.index = heat_front.index.map(region_display_names)
heat_front.columns = month_names

front_start_year = int(front_tbl["front_start"].dt.year.min())
front_end_year = int(front_tbl["front_start"].dt.year.max())

fig, ax = plt.subplots(figsize=(4, 1.5))
sns.heatmap(heat_front,ax=ax,cmap="gray_r",annot=True,fmt=".0f",linewidths=0.5,
    linecolor="white",cbar=False,
)
ax.set_title(f"Frontal passages by month ({front_start_year}-{front_end_year})")
ax.set_xlabel("")
ax.set_ylabel("")
fig.tight_layout()

finished_plot(fig,"../plots/manuscript/fig1_seasonal_fronts_heatmap.jpg")



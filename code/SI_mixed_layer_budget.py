#%%
"""Mixed-layer heat-flux budget diagnostic.

This script estimates the heat flux needed to generate a given temperature
change over a specified time window for different mixed-layer depths (MLD).
The calculations are based on the simple bulk mixed-layer relationship:

    Q = rho * cp * h * dT / tau

where rho is seawater density, cp is specific heat capacity, h is mixed-layer
thickness, dT is the temperature change, and tau is the time interval.
"""

import os
import sys
from matplotlib.lines import Line2D
# Startup config for the project environment.
STARTUP_2026_coldfronts_mhws = "../utils/2026_coldfronts_mhws_startup.py"
exec(open(STARTUP_2026_coldfronts_mhws).read())

sys.path.append("../")
from utils.plot_utils import finished_plot


# ---------------------------------------------------------------------------
# Simple mixed-layer heat budget estimate from a single MLD case
# ---------------------------------------------------------------------------

mld = 50  # mixed-layer depth (m)

dTs = np.arange(0, 1.1, 0.05)  # temperature change range (degC)
q = np.array([1025 * 3990 * mld * dT for dT in dTs])
# Total heat required to warm a 50-m mixed layer by dT, in J/m^2.
dt = np.arange(6, 5 * 24, 1)  # time window in hours

# Convert heat to equivalent surface heat flux in W/m^2.
F = [q / (t * 3600) for t in dt]

fig, ax = plt.subplots(figsize=(5, 3))
cc = plt.pcolormesh(
    dTs,
    dt / 24,
    F,
    vmin=0,
    vmax=1000,
    cmap=plt.cm.get_cmap("Spectral_r", 30),
)
# plt.contour(np.arange(0, 1.1, 0.1), dt / 24, F, colors='k')
plt.colorbar(cc)
ax.set_ylabel("Time window (days)")
ax.set_xlabel("Temperature change (°C)")


# ---------------------------------------------------------------------------
# Bulk mixed-layer heat-flux maps for several MLD values
# ---------------------------------------------------------------------------

rho = 1025.0  # kg/m^3
cp = 3990.0   # J/(kg K)

# Temperature changes and timescales plotted in the diagnostic figures.
dTs = np.arange(0, 1.05, 0.05)       # degC
taus = np.arange(6, 5 * 24 + 1, 1) / 24  # days
mlds = [10, 20, 50, 100]             # mixed-layer depth (m)

# Convert the time axis to seconds for the heat-flux calculation.
tau_sec = taus * 86400.0

# Grid of dT and tau values for contour plots.
D, T = np.meshgrid(dTs, tau_sec)  # D = delta T, T = time in seconds

fig, axes = plt.subplots(2, 2, figsize=(8, 6), sharex=True, sharey=True)

vmin, vmax = 0, 1500
levels = np.linspace(vmin, vmax, 31)

for ax, h in zip(axes.flat, mlds):
    # Required heat flux to produce a given dT over the time tau.
    F = rho * cp * h * D / T  # W/m^2

    cf = ax.contourf(dTs, taus, F, levels=levels, cmap="Spectral_r", extend="max")
    ax.set_title(f"MLD = {h} m")
    ax.set_xlabel(r"$\Delta T$ ($^\circ$C)")
    ax.set_ylabel(r"$\tau$ (days)")
    ax.grid(alpha=0.2)

cbar = plt.colorbar(cf, ax=axes[:], shrink=0.9)
cbar.set_label(r"Required heat flux (W m$^{-2}$)")

plt.tight_layout()
plt.show()


# ---------------------------------------------------------------------------
# Focused comparison figure for a representative MLD range - for manuscript
# ---------------------------------------------------------------------------

font_for_pres()
fig, ax = plt.subplots(figsize=(6, 4))

# Filled contours for MLD = 20 m.
F20 = rho * cp * 20 * D / T
cf = ax.contourf(
    dTs,
    taus,
    F20,
    levels=np.arange(0, 1600, 100),
    cmap="Spectral_r",
    extend="max",
)
cs = ax.contour(
    dTs,
    taus,
    F20,
    levels=np.arange(0, 1700, 400),
    colors='k',
    alpha=0.3,
)
ax.clabel(cs, fmt='%d')
plt.colorbar(cf, label=r"Heat flux (W m$^{-2}$)", shrink=0.6)

# Overplot black contours for MLD = 100 m for reference.
F100 = rho * cp * 100 * D / T
cs = ax.contour(
    dTs,
    taus,
    F100,
    levels=np.arange(0, 1600, 400),
    colors='k',
    linewidths=1,
    linestyles='--',
)
ax.clabel(cs, fmt='%d')

ax.set_xlabel(r"$\Delta T$ ($^\circ$C)")
ax.set_ylabel(r"$\tau$ (days)")

legend_lines = [
    Line2D(
        [0], [0],
        color='k',
        lw=1.5,
        linestyle='-',
        alpha=0.4,
        label='MLD = 20 m',
    ),
    Line2D(
        [0], [0],
        color='k',
        lw=1.5,
        linestyle='--',
        label='MLD = 100 m',
    ),
]

ax.legend(handles=legend_lines, loc='upper right', frameon=True)


# save figure
finished_plot(fig,'../plots/manuscript/SI_mixed_layer_budget_estimate.png')

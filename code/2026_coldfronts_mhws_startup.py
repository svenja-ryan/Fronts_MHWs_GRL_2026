import os
#standart libraries:
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator
# Cartopy for geospatial plotting
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy
import pandas as pd
import string
import sys
import xarray as xr
import dask
import scipy.io as sc
import pickle
import warnings
warnings.filterwarnings("ignore")


# names for locations
names = ['JordanBasin','ScotianShelf1','ScotianShelf2','Pioneer','SlopeN','SlopeS','MAB_North','MAB_south']
# lat, lon for location
# closest ERA 5 locations to initial locations --> provided by Rhys
lat_fronts = [43.5,43,44.5,40.25,41.5,39,40.75,38.5]
lon_fronts = np.array([292,295,296.5,289,295,289,288,285.5])-360

# Dictionary keyed by location name for convenient coordinate lookup
front_loc = {
    name: {"lat": lat, "lon": lon}
    for name, lat, lon in zip(names, lat_fronts, lon_fronts)
}


# convert matlab time to datetime
def matlab2datetime(matlab_datenum):
    day = dt.datetime.fromordinal(int(matlab_datenum))
    dayfrac = dt.timedelta(days=matlab_datenum%1) - dt.timedelta(days = 366)
    return day + dayfrac

#----------------------------------------------------------------------
# map projections
#----------------------------------------------------------------------
proj_lamb = ccrs.LambertConformal(
    central_longitude=-65,
    central_latitude=45,
    standard_parallels=(35, 55),
)

# ---------------------------------------------------------------------
# Load Bathymetry Data
# ---------------------------------------------------------------------
# Subset East Coast: longitude –85 to –40, latitude 25 to 55
bathy = (
    xr.open_dataset('/vast/clidex/data/bathymetry/ETOPO1/ETOPO1_Bed_g_gmt4.grd')
      .sel(x=slice(-85, -30), y=slice(25, 55))
      .load()
)


# ---------------------------------------------------------------------
# Load pioneer mooring bulk stratification MAT file into an xarray Dataset
# ---------------------------------------------------------------------
def load_pioneer_mooring(filepath =  "/srv/data/Pioneer_NE/Data_BulkStratInshore.mat"): # -> xr.Dataset:
    """Load the Pioneer mooring bulk stratification MAT file into an xarray Dataset."""
    mat = scipy.io.loadmat(str(filepath), squeeze_me=True, struct_as_record=False)
    strat = mat["strat_bulk"]

    time = [matlab2datetime(t) for t in strat.time]
    ddens = np.asarray(strat.ddens, dtype="float64")
    pres = np.asarray(strat.pres, dtype="float64")

    ds = xr.Dataset(
        data_vars={
            "ddens": (("time",), ddens),
            "pres": (("pressure",), pres),
        },
        coords={
            "time": time,
            "pressure": pres,
        },
        attrs={
            "source": str(filepath),
            "description": "Bulk stratification from Pioneer mooring MAT file",
        },
    )
    return ds



# ---------------------------------------------------------------------
# plot map
# ---------------------------------------------------------------------
def plot_map_NWA(ax,extent,c=np.arange(0,6000,500),plotbathy=None,gulfstream=False):
    """
    Produces map with coastlines using cartopy.
    
    INPUT:
    ax       : axis handle 
    extent   : [lon1,lon2,lat1,lat2]
    
    OUTPUT:
    gl       : grid handle
    
    """
    
    ## load bathymetry
    bathy = xr.open_dataset('/vast/clidex/data/bathymetry/ETOPO2v2/ETOPO2v2c_f4.nc').sel(x=slice(-85,-40),y=slice(25,55)).load()

    
    # projection
    proj = ccrs.PlateCarree()
    
    ## plot bathymetry
    if plotbathy=='contourf':
        ax.contourf(bathy.x,bathy.y,(bathy.z*(-1)),levels=c,cmap=plt.get_cmap('Blues',len(c)),transform=ccrs.PlateCarree())
    elif plotbathy=='contour':
        ax.contour(bathy.x,bathy.y,(bathy.z*(-1)),levels=c,colors='gray',transform=ccrs.PlateCarree(),linewidths=1)
        
    ## add mean Gulf Stream path
    if gulfstream:
        gs = sc.loadmat('/srv/data/GS_monthly_yearly_two_sat_contours_1993_to_2018.mat')
        gs_time = gs['time_monthly']
        gs_lon = gs['lon_cell'][0,0][0,:]
        gs_lat = gs['lat_cell'][0,0][0,:]
        plt.plot(gs_lon,gs_lat,transform=ccrs.PlateCarree(),color='k',linewidth=3,alpha=0.3)

    ## formatting
    ax.set_extent(extent)
#     ax.coastlines(resolution='50m',color='gray')
    gl = ax.gridlines(crs=proj,draw_labels=True,linestyle='-',alpha=0.3,
                      xlocs=range(-120,40,5),ylocs=range(0,80,5),x_inline=False)
    ax.add_feature(cartopy.feature.GSHHSFeature(scale='intermediate'), facecolor='lightgray',edgecolor='gray')
#     gl.xlocator = mticker.FixedLocator(np.arange(-80,-40,5))
#     gl.ylocator = mticker.FixedLocator([10, 20, 30, 40, 50])
    gl.rotate_labels = False
    gl.top_labels = False
    gl.right_labels = False
    
    return gl


#
#------------------------------------------------------------------------
# If a figure name is defined, save the figure to that file. Otherwise, display the figure on screen.
def finished_plot (fig, fig_name=None, dpi=300):

    if fig_name is not None:
        print('Saving ' + fig_name)
        fig.savefig(fig_name, dpi=dpi,bbox_inches='tight')
    else:
        fig.show()


# my own libraries:
#import m_general as M

#import AA_plot_base as AA
def json_load(name, path, verbose=False):
    import json
    full_name= (os.path.join(path,name+ '.json'))

    with open(full_name, 'r') as ifile:
        data=json.load(ifile)
    if verbose:
        print('loaded from: ',full_name)
    return data

# mconfig=json_load('config','/home/mhell/2021_ICESat2_tracks/config/')

# # add project depenent libraries
# sys.path.append(mconfig['paths']['local_script'])
# sys.path.append(mconfig['paths']['local_script'] +'/ICEsat2_SI_tools/')

# import m_colormanager_ph3 as M_color
# import m_tools_ph3 as MT
# import m_general_ph3 as M

#load colorscheme
# col=M_color.color(path=mconfig['paths']['analysis']+'../config/', name='color_def')


lstrings =iter([i+') ' for i in list(string.ascii_lowercase)])
# define journal fig sizes
# fig_sizes = mconfig['fig_sizes']['AMS']


SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12
#csfont = {'fontname':'Comic Sans MS'}
legend_properties = {'weight':'bold'}
#font.family: sans-serif
#font.sans-serif: Helvetica Neue

#import matplotlib.font_manager as font_manager
#font_dirs = ['/home/mhell/HelveticaNeue/', ]
#font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
#font_list = font_manager.createFontList(font_files)
#font_manager.fontManager.ttflist.extend(font_list)

plt.rc('font', size=SMALL_SIZE, serif='Helvetica Neue', weight='normal')          # controls default text sizes
#plt.rc('font', size=SMALL_SIZE, serif='DejaVu Sans', weight='light')
plt.rc('text', usetex='false')
plt.rc('axes', titlesize=MEDIUM_SIZE, labelweight='normal')     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE, labelweight='normal') #, family='bold')    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE, frameon=False)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE, titleweight='bold', autolayout=True) #, family='bold')  # fontsize of the figure title

#figure.autolayout : False
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42


plt.rc('path', simplify=True)

plt.rcParams['figure.figsize'] = (10, 8)#(20.0, 10.0) #inline


### TICKS
# see http://matplotlib.org/api/axis_api.html#matplotlib.axis.Tick
#xtick.top            : False   # draw ticks on the top side
#xtick.bottom         : True   # draw ticks on the bottom side
#xtick.major.size     : 3.5      # major tick size in points
#xtick.minor.size     : 2      # minor tick size in points
#xtick.major.width    : .8    # major tick width in points
#xtick.minor.width    : 0.6    # minor tick width in points
#xtick.major.pad      : 3.5      # distance to major tick label in points
#xtick.minor.pad      : 3.4      # distance to the minor tick label in points
#xtick.color          : k      # color of the tick labels
#xtick.labelsize      : medium # fontsize of the tick labels
#xtick.direction      : out    # direction: in, out, or inout
#xtick.minor.visible  : False  # visibility of minor ticks on x-axis
#xtick.major.top      : True   # draw x axis top major ticks
#xtick.major.bottom   : True   # draw x axis bottom major ticks
#xtick.minor.top      : True   # draw x axis top minor ticks
#xtick.minor.bottom   : True   # draw x axis bottom minor ticks

#ytick.left           : True   # draw ticks on the left side
#ytick.right          : False  # draw ticks on the right side
#ytick.major.size     : 3.5      # major tick size in points
#ytick.minor.size     : 2      # minor tick size in points
#ytick.major.width    : 0.8    # major tick width in points
#ytick.minor.width    : 0.6    # minor tick width in points
#ytick.major.pad      : 3.5      # distance to major tick label in points
#ytick.minor.pad      : 3.4      # distance to the minor tick label in points
#ytick.color          : k      # color of the tick labels
#ytick.labelsize      : medium # fontsize of the tick labels
#ytick.direction      : out    # direction: in, out, or inout
#ytick.minor.visible  : False  # visibility of minor ticks on y-axis
#ytick.major.left     : True   # draw y axis left major ticks
#ytick.major.right    : True   # draw y axis right major ticks
#ytick.minor.left     : True   # draw y axis left minor ticks
#ytick.minor.right    : True   # draw y axis right minor ticks


plt.rc('xtick.major', size= 4, width=1 )
plt.rc('ytick.major', size= 3.8, width=1 )

#axes.facecolor      : white   # axes background color
#axes.edgecolor      : black   # axes edge color
#axes.linewidth      : 0.8     # edge linewidth
#axes.grid           : False   # display grid or not
#axes.titlesize      : large   # fontsize of the axes title
#axes.titlepad       : 6.0     # pad between axes and title in points
#axes.labelsize      : medium  # fontsize of the x any y labels
#axes.labelpad       : 4.0     # space between label and axis
#axes.labelweight    : normal  # weight of the x and y labels
#axes.labelcolor     : black

plt.rc('axes', labelsize= MEDIUM_SIZE, labelweight='normal')




# axes.spines.left   : True   # display axis spines
# axes.spines.bottom : True
# axes.spines.top    : True
# axes.spines.right  : True
plt.rc('axes.spines', top= False, right=False )



def font_for_print():

    SMALL_SIZE = 6
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 10
    #csfont = {'fontname':'Comic Sans MS'}
    legend_properties = {'weight':'bold'}
    #font.family: sans-serif
    #font.sans-serif: Helvetica Neue

    #import matplotlib.font_manager as font_manager
    #font_dirs = ['/home/mhell/HelveticaNeue/', ]
    #font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
    #font_list = font_manager.createFontList(font_files)
    #font_manager.fontManager.ttflist.extend(font_list)

    plt.rc('font', size=SMALL_SIZE, serif='Helvetica Neue', weight='normal')          # controls default text sizes
    #plt.rc('font', size=SMALL_SIZE, serif='DejaVu Sans', weight='light')
    plt.rc('text', usetex='false')
    plt.rc('axes', titlesize=MEDIUM_SIZE, labelweight='normal')     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE, labelweight='normal') #, family='bold')    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE, frameon=False)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE, titleweight='bold', autolayout=True) #, family='bold')  # fontsize of the figure title

    #figure.autolayout : False
    #matplotlib.rcParams['pdf.fonttype'] = 42
    #matplotlib.rcParams['ps.fonttype'] = 42


    #plt.rc('xtick.major', size= 4, width=1 )
    #plt.rc('ytick.major', size= 3.8, width=1 )


    plt.rc('axes', labelsize= SMALL_SIZE, labelweight='normal')

def font_for_pres():

    SMALL_SIZE = 10
    MEDIUM_SIZE = 12
    BIGGER_SIZE = 14
    #csfont = {'fontname':'Comic Sans MS'}
    legend_properties = {'weight':'bold'}
    #font.family: sans-serif
    #font.sans-serif: Helvetica Neue

    #import matplotlib.font_manager as font_manager
    #font_dirs = ['/home/mhell/HelveticaNeue/', ]
    #font_files = font_manager.findSystemFonts(fontpaths=font_dirs)
    #font_list = font_manager.createFontList(font_files)
    #font_manager.fontManager.ttflist.extend(font_list)

    plt.rc('font', size=SMALL_SIZE, serif='Helvetica Neue', weight='normal')          # controls default text sizes
    #plt.rc('font', size=SMALL_SIZE, serif='DejaVu Sans', weight='light')
    plt.rc('text', usetex='false')
    plt.rc('axes', titlesize=MEDIUM_SIZE, labelweight='normal')     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE, labelweight='normal') #, family='bold')    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE, frameon=False)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE, titleweight='bold', autolayout=True) #, family='bold')  # fontsize of the figure title

    #figure.autolayout : False
    #matplotlib.rcParams['pdf.fonttype'] = 42
    #matplotlib.rcParams['ps.fonttype'] = 42


    #plt.rc('xtick.major', size= 4, width=1 )
    #plt.rc('ytick.major', size= 3.8, width=1 )


    plt.rc('axes', labelsize= SMALL_SIZE, labelweight='normal')



# add project depenent libraries
#sys.path.append(config['paths']['local_script'])


def map_init(ax=None,extent=[-76, -63, 36, 45]):
    """
    Plot one or more xarray datasets (`datasets`) of a scalar variable `var`
    over the Mid‐Atlantic Bight region. Adds GSHHS coastline, 100 m bathymetry,
    gridlines, and an optional colorbar if any data were plotted.
    """
    font_for_pres()
    # Create new figure & PlateCarree axis if needed
    if ax is None:
        fig, ax = plt.subplots(
            figsize=(6, 5),
            subplot_kw=dict(projection=ccrs.PlateCarree())
        )
    else:
        fig = ax.figure  # <-- IMPORTANT
        
    if extent is None:
        ax.set_extent(extent,crs=ccrs.PlateCarree())
    else:
        ax.set_extent(extent,crs=ccrs.PlateCarree())

    # Track last mappable (for colorbar)
    cc = None

    # Add gridlines
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(), draw_labels=True,
        linewidth=1, color='lightgray', alpha=0.5, linestyle='--'
    )
    gl.top_labels = False
    gl.right_labels = False

    # Add GSHHS coastline
    ax.add_feature(
        cfeature.GSHHSFeature(scale='high'),
        facecolor='lightgray', edgecolor='None'
    )

    # Contour 100 m bathymetry
    ax.contour(
        bathy.x, bathy.y, bathy.z * (-1),
        levels=[50,75,100,200, 1000], colors='k',
        transform=ccrs.PlateCarree(),
        linewidths=1, alpha=0.3
    )

    return fig, ax



import os
import pickle

#
#------------------------------------------------------------------------
# 1) save pickle file

def pickle_save(name, path, data, verbose=True):
    if not os.path.exists(path):
        os.makedirs(path)
    full_name= (os.path.join(path,name+ '.npy'))


    with open(full_name, 'wb') as f2:
        pickle.dump(data, f2)
    if verbose:
        print('save at: ',full_name)
        
#
#------------------------------------------------------------------------
# 1) load pickle file  

def pickle_load(name, path, verbose=True):  
    #if not os.path.exists(path):
    #    os.makedirs(path)
    full_name= (os.path.join(path,name+ '.npy'))

    with open(full_name, 'rb') as f:
        data=pickle.load(f)

    if verbose:
        print('load from: ',full_name)
    return data
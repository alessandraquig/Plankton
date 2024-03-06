import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LogNorm, CenteredNorm
from matplotlib.ticker import FuncFormatter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import os
import numpy as np
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
matplotlib.use('Agg')

from matplotlib import rc
this_rc_params = {
    "text.usetex": True,
    "font.family": "roman"
}
plt.rcParams.update(this_rc_params)
#matplotlib.rcParams['figure.dpi'] = 400

def plot_SiNO3(layer):
    if layer == 'depth':
        path = "Data/nc/n_si_ratio_200m.nc"
    elif layer == 'surf':
        path = "Data/nc/n_si_ratio_surf.nc"

    # Read the netCDF file
    data = nc.Dataset(path)
    var_name = os.path.basename(path).split(".")[0]  # Check variable name first

    # Extract the latitude and longitude variables
    lat_var = data.variables.get('y')
    lat = lat_var[:]
    lon_var = data.variables.get('x') or data.variables.get('longitude')
    lon = lon_var[:]

    # Extract the data variable you want to plot
    var = data.variables[var_name][:]
    print(var.min(), var.max())

    # Create a figure and axes with a specific projection
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

    # Plot the data on a latitude and longitude scale
    im = ax.pcolormesh(lon, lat, var, transform=ccrs.PlateCarree(), cmap='PiYG', norm=CenteredNorm(vcenter=1))

    # Set the extent of the map to match your data
    ax.set_extent([-160, -70, -60, 0], crs=ccrs.PlateCarree())

    # Add parallels and meridians
    ax.gridlines(draw_labels=False, linewidth=0.5, color='grey', alpha=0.5, linestyle='-')
    ax.set_xticks(np.arange(-160, -60, 10), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(-60, 10, 10), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())

    # Mask out land
    ax.add_feature(cfeature.LAND, zorder=1, facecolor='w')

    # Add coastlines
    ax.coastlines('50m')
    #ax.add_feature(cfeature.COASTLINE, zorder=2)

    # Add a colorbar
    cbar = plt.colorbar(im, ax=ax)

    # Set the title and labels
    ax.set_title(r'Si/NO$^{3}$ Ratio')
    ax.set_xlabel(r'Longitude')
    ax.set_ylabel(r'Latitude')

    return im, ax, cbar

if __name__ == "__main__":
    layer = 'surf'

    im, ax, cbar = plot_SiNO3(layer)

    if layer == 'depth':
        plt.savefig('Output/SiNO3_meso.pdf', format='pdf')
    elif layer == 'surf':
        plt.savefig('Output/SiNO3_epi.pdf', format='pdf')


    # Close the plot
    plt.close()
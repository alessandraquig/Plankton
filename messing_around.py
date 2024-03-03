import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import numpy as np
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from matplotlib.colors import LogNorm

matplotlib.use('Agg')

from matplotlib import rc
this_rc_params = {
    "text.usetex": True,
    "font.family": "roman"
}
plt.rcParams.update(this_rc_params)

def plot_chlorophyll():
    path = "Data/chlorophyll_mean.nc"
    data = nc.Dataset(path)
    var_name = "chlorophyll"

    # Extract the latitude and longitude variables
    lat_var = data.variables.get('latitude')
    lat = lat_var[:]
    lon_var = data.variables.get('longitude')
    lon = lon_var[:]

    # Extract the data variable you want to plot
    var = data.variables[var_name][:]
    print(f'min = {var.min()}, max = {var.max()}, mean = {var.mean()}')

    # Create a figure and axes with a specific projection
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

    # Plot the data on a latitude and longitude scale using the custom colormap
    im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap='Greens', norm=LogNorm(), levels=6)

    # Set the extent of the map to match your data
    ax.set_extent([-180, -65, -70, 0], crs=ccrs.PlateCarree())

    # Add parallels and meridians
    ax.gridlines(draw_labels=False, linewidth=0.5, color='grey', alpha=0.5, linestyle='-')
    ax.set_xticks(np.arange(-180, -55, 10), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(-70, 10, 10), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())

    # Mask out land
    ax.add_feature(cfeature.LAND, zorder=1, facecolor='w')

    # Add coastlines
    ax.coastlines('50m')

    # Add a colorbar
    cbar = plt.colorbar(im, ax=ax, extend='max', ticks=np.linspace(var.min(), var.max(), 6))

    # Set the title and labels
    ax.set_title(r'Surface Chlorophyll Concentration ($mg/m^{3}$)')
    ax.set_xlabel(r'Longitude')
    ax.set_ylabel(r'Latitude')

    return im, ax, cbar

if __name__ == "__main__":

    im, ax, cbar = plot_chlorophyll()
    
    # Save the plot as a tif file
    plt.savefig('Output/chlorophyll_mean.tif', format='tif')

    # Close the plot
    plt.close()

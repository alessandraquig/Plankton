import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import os
import numpy as np
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
matplotlib.use('qtagg')

from matplotlib import rc
this_rc_params = {
    "text.usetex": True,
    "font.family": "roman"
}
plt.rcParams.update(this_rc_params)
#matplotlib.rcParams['figure.dpi'] = 400

def plot_SST():

    # Read the netCDF file
    path = "Data/SST_mean.nc"
    data = nc.Dataset(path)
    var_name = os.path.basename(path).split(".")[0]
    if var_name == "SST_mean":
        var_name = "temperature"
    print(f"Variable name: {var_name}")

    # Extract the latitude and longitude variables
    lat_var = data.variables.get('lat') or data.variables.get('latitude')
    if lat_var is not None:
        lat = lat_var[:]
    else:
        raise ValueError("Latitude variable not found in the netCDF file.")
    lon_var = data.variables.get('lon') or data.variables.get('longitude')
    if lon_var is not None:
        lon = lon_var[:]
    else:
        raise ValueError("Longitude variable not found in the netCDF file.")

    # Extract the data variable you want to plot
    var = data.variables[var_name][:]
    print(var.min(), var.max())

    # Create a figure and axes with a specific projection
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

    # Creating levels
    #norm = BoundaryNorm(levels=np.linspace(), ncolors=cmap.N, clip=True)

    # Plot the data on a latitude and longitude scale
    #im = ax.pcolormesh(lon, lat, var, transform=ccrs.PlateCarree(),  cmap='coolwarm')

    # Plot the data on a latitude and longitude scale
    im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap='coolwarm')

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
    #ax.add_feature(cfeature.COASTLINE, zorder=2)

    # Add a colorbar
    cbar = plt.colorbar(im, ax=ax)

    # Set the title and labels
    ax.set_title(r'Sea Surface Temperature ($^{\circ}$C)')
    ax.set_xlabel(r'Longitude')
    ax.set_ylabel(r'Latitude')

    return im, ax, cbar

if __name__ == '__main__':

    im, ax, cbar = plot_SST()

    # Save the plot as a tif file
    plt.savefig('Output/SST.tif', format='tif')

    # Close the plot
    plt.close()
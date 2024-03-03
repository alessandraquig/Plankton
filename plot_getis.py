import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
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

def plot_getis(plankton, layer):
    """
    :param plankton: "phyto" or "zoo"
    :param layer: "surf" or "depth"
    :return: fig, ax, cbar
    """
    path = f"Data/{plankton}getis{layer}.nc"

    data = nc.Dataset(path)
    var_name = os.path.basename(path).split(".")[0] # Check this

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
    var[var <= 0] = 1e-2
    print(var.min(), var.max())


    # Create a figure and axes with a specific projection
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

    # Norm the data so 1 is neutral
    #norm=colors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=3.5)
    norm = matplotlib.colors.LogNorm()

    # Plot the data on a latitude and longitude scale
    im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap='PiYG', norm=norm)#, norm=matplotlib.colors.LogNorm(), levels=np.logspace(np.log10(var.min()), np.log10(var.max()), 10), extend='max')

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
    if layer == "surf":
        layer_name = "Epipelagic"
    if layer == "depth":
        layer_name = "Mesopelagic"
    ax.set_title(rf'Getis Ord')
    #ax.set_xlabel(r'Longitude')
    #ax.set_ylabel(r'Latitude')

    return im, ax, cbar

if __name__ == "__main__":

    layer = 'surf'
    plankton = 'phyto'
    if layer == "surf":
        layer_name = "epi"
    if layer == "depth":
        layer_name = "meso"

    im, ax, cbar = plot_getis('phyto', 'surf')
    ax.set_title(rf'{layer_name.capitalize()}pelagic {plankton.capitalize()}plankton Getis Ord')

    # Save the plot as a tif file
    plt.savefig(f'Output/{layer_name}_{plankton}_getis.tif', format='tif')

    # Close the plot
    plt.close()
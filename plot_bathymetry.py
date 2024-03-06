import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import numpy as np
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

matplotlib.use('Agg')

from matplotlib import rc
this_rc_params = {
    "text.usetex": True,
    "font.family": "roman"
}
plt.rcParams.update(this_rc_params)

def plot_bathymetry():

    # Read the netCDF file
    path = "Data/bathymetry.nc"
    data = nc.Dataset(path)
    var_name = 'elevation'

    # Extract the latitude and longitude variables and downsample by a factor of 10
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

    # Downsample the data
    downsample_factor = 10  # Adjust this value as needed
    downsampled_lon = lon[::downsample_factor]
    downsampled_lat = lat[::downsample_factor]
    downsampled_var = var[::downsample_factor, ::downsample_factor]
    print("Passed downsampling")

    # Create a figure and axes with a specific projection
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

    # Creating levels
    #norm = BoundaryNorm(levels=np.linspace(), ncolors=cmap.N, clip=True)

    # Plot the data on a latitude and longitude scale
    #im = ax.pcolormesh(lon, lat, var, transform=ccrs.PlateCarree(),  cmap='coolwarm')

    # Plot the data on a latitude and longitude scale
    im = ax.pcolormesh(downsampled_lon, downsampled_lat, downsampled_var, transform=ccrs.PlateCarree(), cmap='viridis', vmax=0)
    print('Created image')

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
    ax.set_title('Depth ($m$)')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    return im, ax, cbar

if __name__ == '__main__':

    im, ax, cbar = plot_bathymetry()
    print('Plotted bathymetry')

    # Save the plot as a tif file
    plt.savefig('Output/bathymetry_map.pdf', format='pdf')
    print('Saved figure')

    # Close the plot
    plt.close()
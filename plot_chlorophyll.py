import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import numpy as np
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from matplotlib.colors import LinearSegmentedColormap

matplotlib.use('qtagg')

from matplotlib import rc
this_rc_params = {
    "text.usetex": True,
    "font.family": "roman"
}
plt.rcParams.update(this_rc_params)

# Read the netCDF file
path = "Data/chlorophyll_mean.nc"
data = nc.Dataset(path)
var_name = "chlorophyll"
print(f"Variable name: {var_name}")

# Extract the latitude and longitude variables
lat_var = data.variables.get('latitude')
lat = lat_var[:]
lon_var = data.variables.get('longitude')
lon = lon_var[:]

# Extract the data variable you want to plot
var = data.variables[var_name][:]
print(f'min = {var.min()}, max = {var.max()}, mean = {var.mean()}')

# Define custom levels for contourf
levels = [0, 0.1, 0.3, 0.6, 1, 2, 5]

# Create a custom colormap emphasizing low values
colors = ['#f7fcf5', '#a1d99b', '#74c476', '#41ab5d', '#238b45', '#127031', '#00441b']
cmap = LinearSegmentedColormap.from_list("custom_green", colors)

# Create a figure and axes with a specific projection
fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

# Plot the data on a latitude and longitude scale using the custom colormap
im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap=cmap, levels=levels) #, extend='max')

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
cbar = plt.colorbar(im, ax=ax)

# Set the title and labels
ax.set_title(r'Surface Chlorophyll Concentration ($mg/m^{3}$)')
ax.set_xlabel(r'Longitude')
ax.set_ylabel(r'Latitude')

# Show the plot
plt.show()

# Save the plot as a tif file
plt.savefig('Output/chlorophyll_mean.tif', format='tif')

# Close the plot
plt.close()
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

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

# Create a figure and axes with a specific projection
fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

# Creating levels
boundaries = [0, 0.5, 1, 2, 3, 4, 5]

# Plot the data on a latitude and longitude scale
im = ax.pcolormesh(lon, lat, var, transform=ccrs.PlateCarree(), cmap='Greens', norm=BoundaryNorm(boundaries, len(boundaries)-1))

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
ax.set_title('Surface Chlorophyll Concentration (mg/m^3)')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Show the plot
plt.show()

# Save the plot as a tif file
plt.savefig('Output/chlorophyll_mean.tif', format='tif')

# Close the plot
plt.close()
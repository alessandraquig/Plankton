import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import os
import numpy as np

# Read the netCDF file
path = "Data/phytonestednesssurf.nc"
data = nc.Dataset(path)
var_name = os.path.basename(path).split(".")[0]
print(f"Variable name: {var_name}")

# Extract the latitude and longitude variables
lat = data.variables['lat'][:]
lon = data.variables['lon'][:]

# Extract the data variable you want to plot
var = data.variables[var_name][:]

# Create a figure and axes with a specific projection
fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

# Plot the data on a latitude and longitude scale
im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), levels=np.linspace(0, 1, 10 + 1),  cmap='viridis')

# Mask out land
ax.add_feature(cfeature.LAND, zorder=1, facecolor='w')

# Add coastlines
ax.coastlines('10m')
#ax.add_feature(cfeature.COASTLINE, zorder=2)

# Add a colorbar
cbar = plt.colorbar(im, ax=ax)

# Set the title and labels
ax.set_title('Phytoplankton Surface Nestedness')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Show the plot
plt.show()
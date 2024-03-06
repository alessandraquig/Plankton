import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import os
import numpy as np
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Read the netCDF file
path = "Data/nc/n_si_ratio_surf.nc"
data = nc.Dataset(path)
var_name = os.path.basename(path).split(".")[0]  # Check variable name first

# Extract the latitude and longitude variables
lat_var = data.variables.get('y')
lat = lat_var[:]
lon_var = data.variables.get('x') or data.variables.get('longitude')
lon = lon_var[:]

# Extract the data variable you want to plot
var = data.variables[var_name][:]

# Check for NaNs in the data
nan_mask = np.isnan(var)

# Replace NaNs with 1000
var[nan_mask] = 1000

# Create a figure and axes with a specific projection
fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

# Define the MidPointLogNorm class
class MidPointLogNorm(LogNorm):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        LogNorm.__init__(self,vmin=vmin, vmax=vmax, clip=clip)
        self.midpoint=midpoint
    def __call__(self, value, clip=None):
        value = np.ma.masked_invalid(value)
        x, y = [np.log(self.vmin), np.log(self.midpoint), np.log(self.vmax)], [0, 0.5, 1]
        print(f'Min value = {value.min()}, max value = {value.max()}')
        z = np.ma.masked_array(np.interp(np.log(value), x, y))
        print(z.min(), z.max())
        return z

# Plot the data on a latitude and longitude scale
im = ax.pcolormesh(lon, lat, var, transform=ccrs.PlateCarree(), cmap='PiYG', norm=MidPointLogNorm(vmin=0.0001, vmax=10, midpoint=1))

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

# Add a colorbar
cbar = plt.colorbar(im, ax=ax)

# Set the title and labels
ax.set_title(r'Si/NO$^{3}$ Ratio')
ax.set_xlabel(r'Longitude')
ax.set_ylabel(r'Latitude')

# Save the plot as a tif file
plt.savefig('Output/SiNO3.tif', format='tif')

# Close the plot
plt.close()

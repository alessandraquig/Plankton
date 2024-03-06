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
print(var.min(), var.max())

# Check for NaNs in the data
nan_mask = np.isnan(var)

# Replace NaNs with 1000
var[nan_mask] = 1000

#Taking log of the ratio
var = np.log10(var)


# Create a figure and axes with a specific projection
fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

# Creating levels
#norm = BoundaryNorm(levels=np.linspace(), ncolors=cmap.N, clip=True)

# Plot the data on a latitude and longitude scale
#im = ax.pcolormesh(lon, lat, var, transform=ccrs.PlateCarree(),  cmap='coolwarm')



# Plot the data on a latitude and longitude scale
im = ax.pcolormesh(lon, lat, var, transform=ccrs.PlateCarree(), cmap='PiYG', norm=CenteredNorm(vcenter=0))

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

# Define your custom tick formatter function to undo the log transformation
def custom_formatter(x, pos):
    return f'{10**x:.0f}'  # Undo the log transformation

# Create your lambda function
custom_formatter_lambda = lambda x, pos: f'{10**x:.0f}'

# Add a colorbar
cbar = plt.colorbar(im, ax=ax)
cbar.ax.yaxis.set_major_formatter(FuncFormatter(custom_formatter))
# Or using the lambda function
cbar.ax.yaxis.set_major_formatter(FuncFormatter(custom_formatter_lambda))

# Set the title and labels
ax.set_title(r'Si/NO$^{3}$ Ratio')
ax.set_xlabel(r'Longitude')
ax.set_ylabel(r'Latitude')

# Save the plot as a tif file
plt.savefig('Output/SiNO3.tif', format='tif')

# Close the plot
plt.close()
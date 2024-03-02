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

# Read the netCDF file
layer = "surf"          # "surf" or "depth"
plankton = "phyto"      # "phyto" or "zoo"
path = f"Data/{plankton}nestedness{layer}.nc"


data = nc.Dataset(path)
var_name = os.path.basename(path).split(".")[0] # Check this
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
var[var >= 1] = 1
#var[var < 0] = 1e-8
print(var.min(), var.max())

# Create a figure and axes with a specific projection
fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})

# Creating levels
#norm = BoundaryNorm(levels=np.linspace(), ncolors=cmap.N, clip=True)

# Plot the data on a latitude and longitude scale
#im = ax.pcolormesh(lon, lat, var, transform=ccrs.PlateCarree(),  cmap='coolwarm')

# Plot the data on a latitude and longitude scale
im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap='YlOrRd', levels=np.logspace(-3, 0, 7), norm=matplotlib.colors.LogNorm(vmin=1e-4, vmax=1), extend='min')
#im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap='YlOrRd', norm=matplotlib.colors.LogNorm(vmin=1e-4, vmax=1), extend='min') # levels=np.logspace(np.log10(var.min()), np.log10(var.max()), 12), vmax=1, vmin=1e-4, extend='min')

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
ax.set_title(rf'{layer_name} {plankton.capitalize()}plankton Nestedness')
ax.set_xlabel(r'Longitude')
ax.set_ylabel(r'Latitude')

# Show the plot
plt.show()

# Save the plot as a tif file
plt.savefig(f'Output/{plankton}_nest_{layer}.tif', format='tif')

# Close the plot
plt.close()
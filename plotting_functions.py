import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import numpy as np
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from matplotlib.colors import LogNorm, TwoSlopeNorm
import os

matplotlib.use('Agg')

from matplotlib import rc
this_rc_params = {
    "text.usetex": True,
    "font.family": "roman"
}
plt.rcParams.update(this_rc_params)

def bathymetry(fig, rows, cols, pos=1):

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
    #fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=ccrs.PlateCarree())

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

def DO(layer, fig, rows, cols, pos=1):

    if layer == 'depth':
        path = "Data/nc/DO_200m.nc"
    elif layer == 'surf':
        path = "Data/nc/DO_surf.nc"

    data = nc.Dataset(path)
    var_name = os.path.basename(path).split(".")[0]  # Check variable name first
    print(f"Variable name: {var_name}")

    # Extract the latitude and longitude variables
    lat_var = data.variables.get('y')
    lat = lat_var[:]
    lon_var = data.variables.get('x') or data.variables.get('longitude')
    lon = lon_var[:]

    # Extract the data variable you want to plot
    var = data.variables[var_name][:]
    print(var.min(), var.max())

    # Create a figure and axes with a specific projection
    #fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=ccrs.PlateCarree())

    # Creating levels
    #norm = BoundaryNorm(levels=np.linspace(), ncolors=cmap.N, clip=True)

    # Plot the data on a latitude and longitude scale
    #im = ax.pcolormesh(lon, lat, var, transform=ccrs.PlateCarree(),  cmap='coolwarm')

    # Plot the data on a latitude and longitude scale
    im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap='Reds_r', levels=15) # Reversed colormap to emphasize depletion

    # Set the extent of the map to match your data
    ax.set_extent([-160, -70, -60, 0], crs=ccrs.PlateCarree())

    # Add parallels and meridians
    ax.gridlines(draw_labels=False, linewidth=0.5, color='grey', alpha=0.5, linestyle='-')
    ax.set_xticks(np.arange(-160, -60, 20), crs=ccrs.PlateCarree())
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
    ax.set_title(r'Dissolved Oxygen ($mmol/m^{3}$)')
    ax.set_xlabel(r'Longitude')
    ax.set_ylabel(r'Latitude')

    return im, ax, cbar

def SiNO3(layer, fig, rows, cols, pos=1):
    if layer == 'depth':
        path = "Data/nc/n_si_ratio_200m.nc"
    elif layer == 'surf':
        path = "Data/nc/n_si_ratio_surf.nc"

    # Read the netCDF file
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
    print(np.isnan(np.ma.log(var)))

    # Transform data
    base = 10
    logvar = np.ma.log(var)/np.ma.log(base)

    # Create a figure and axes with a specific projection
    #fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=ccrs.PlateCarree())

    # Plot the data on a latitude and longitude scale
    im = ax.pcolormesh(lon, lat, logvar, transform=ccrs.PlateCarree(), cmap='PiYG', norm=TwoSlopeNorm(vcenter=np.log10(1), vmax=np.log10(10), vmin=np.log10(0.0004)))

    # Set the extent of the map to match your data
    ax.set_extent([-160, -70, -60, 0], crs=ccrs.PlateCarree())

    # Add parallels and meridians
    ax.gridlines(draw_labels=False, linewidth=0.5, color='grey', alpha=0.5, linestyle='-')
    ax.set_xticks(np.arange(-160, -60, 20), crs=ccrs.PlateCarree())
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
    # Place ticks at whole powers of the base.
    cbar.ax.yaxis.set_major_locator(plt.MultipleLocator(1))
    # Define your custom tick formatter function to undo the log transformation
    cbar.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, pos: f'${base:d}^{{{x:.0f}}}$'))

    # Set the title and labels
    ax.set_title(r'Si/NO$^{3}$ Ratio')
    ax.set_xlabel(r'Longitude')
    ax.set_ylabel(r'Latitude')

    return im, ax, cbar

def chlorophyll(fig, rows, cols, pos=1):
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
    #fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=ccrs.PlateCarree())

    # Plot the data on a latitude and longitude scale using the custom colormap
    im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap='Greens', norm=LogNorm(), levels=np.logspace(-2, 1, 7))

    # Set the extent of the map to match your data
    ax.set_extent([-180, -65, -70, 0], crs=ccrs.PlateCarree())

    # Add parallels and meridians
    ax.gridlines(draw_labels=False, linewidth=0.5, color='grey', alpha=0.5, linestyle='-')
    ax.set_xticks(np.arange(-160, -60, 20), crs=ccrs.PlateCarree())
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
    ax.set_title(r'Surface Chlorophyll Concentration ($mg/m^{3}$)')
    ax.set_xlabel(r'Longitude')
    ax.set_ylabel(r'Latitude')

    return im, ax, cbar

def SST(fig, rows, cols, pos=1):

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
    #fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=ccrs.PlateCarree())

    # Creating levels
    #norm = BoundaryNorm(levels=np.linspace(), ncolors=cmap.N, clip=True)

    # Plot the data on a latitude and longitude scale
    #im = ax.pcolormesh(lon, lat, var, transform=ccrs.PlateCarree(),  cmap='coolwarm')

    # Plot the data on a latitude and longitude scale
    im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap='coolwarm')

    # Set the extent of the map to match your data
    ax.set_extent([-160, -70, -60, 0], crs=ccrs.PlateCarree())

    # Add parallels and meridians
    ax.gridlines(draw_labels=False, linewidth=0.5, color='grey', alpha=0.5, linestyle='-')
    ax.set_xticks(np.arange(-160, -60, 20), crs=ccrs.PlateCarree())
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
    ax.set_title(r'Sea Surface Temperature ($^{\circ}$C)')
    ax.set_xlabel(r'Longitude')
    ax.set_ylabel(r'Latitude')

    return im, ax, cbar

def richness(plankton, layer, fig, rows, cols, pos=1):
    """
    :param plankton: "phyto" or "zoo"
    :param layer: "surf" or "depth"
    :return: fig, ax, cbar
    """
    path = f"Data/{plankton}richness{layer}.nc"

    data = nc.Dataset(path)
    var_name = os.path.basename(path).split(".")[0] # Check this
    if plankton == "zoo":
        var_name = "richness"

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
    #fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=ccrs.PlateCarree())

    # Plot the data on a latitude and longitude scale
    im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap='Purples', norm=matplotlib.colors.LogNorm(), levels=np.logspace(np.log10(var.min()), np.log10(var.max()), 10), extend='max')

    # Set the extent of the map to match your data
    ax.set_extent([-160, -70, -60, 0], crs=ccrs.PlateCarree())

    # Add parallels and meridians
    ax.gridlines(draw_labels=False, linewidth=0.5, color='grey', alpha=0.5, linestyle='-')
    ax.set_xticks(np.arange(-160, -60, 20), crs=ccrs.PlateCarree())
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
    cbar.ax.yaxis.set_major_locator(plt.LogLocator())

    # Set the title and labels
    if layer == "surf":
        layer_name = "Epipelagic"
    if layer == "depth":
        layer_name = "Mesopelagic"
    ax.set_title(rf'Richness')
    #ax.set_xlabel(r'Longitude')
    #ax.set_ylabel(r'Latitude')

    return im, ax, cbar

def nestedness(plankton, layer, fig, rows, cols, pos=1):

    path = f"Data/{plankton}nestedness{layer}.nc"


    data = nc.Dataset(path)
    var_name = os.path.basename(path).split(".")[0] # Check this
    if plankton == 'zoo':
        var_name = 'nestedness'
    elif var_name == 'phytonestednessdepth':
        var_name = 'phytonestedness_depth'

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
    #fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=ccrs.PlateCarree())

    # Creating levels
    #norm = BoundaryNorm(levels=np.linspace(), ncolors=cmap.N, clip=True)

    # Plot the data on a latitude and longitude scale
    #im = ax.pcolormesh(lon, lat, var, transform=ccrs.PlateCarree(),  cmap='coolwarm')

    # Plot the data on a latitude and longitude scale
    im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap='YlOrRd', levels=np.logspace(-3, 0, 7), norm=matplotlib.colors.LogNorm(vmin=1e-4, vmax=1), extend='min')
    #im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap='YlOrRd', norm=matplotlib.colors.LogNorm(vmin=1e-4, vmax=1), extend='min') # levels=np.logspace(np.log10(var.min()), np.log10(var.max()), 12), vmax=1, vmin=1e-4, extend='min')

    # Set the extent of the map to match your data
    ax.set_extent([-160, -70, -60, 0], crs=ccrs.PlateCarree())

    # Add parallels and meridians
    ax.gridlines(draw_labels=False, linewidth=0.5, color='grey', alpha=0.5, linestyle='-')
    ax.set_xticks(np.arange(-160, -60, 20), crs=ccrs.PlateCarree())
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
    ax.set_title(rf'Nestedness')

    return im, ax, cbar

def turnover(plankton, layer, fig, rows, cols, pos=1):
    """
    :param plankton: "phyto" or "zoo"
    :param layer: "surf" or "depth"
    :return: fig, ax, cbar
    """
    path = f"Data/{plankton}turnover{layer}.nc"

    data = nc.Dataset(path)
    if plankton == 'zoo':
        var_name = 'turnover'
    elif plankton == 'phyto':
        if layer == 'surf':
            var_name = 'phytturnoversurf'
        elif layer == 'depth':
            var_name = 'phytoturnoverdepth'

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
    var[var > 1] = 1
    print(var.min(), var.max())


    # Create a figure and axes with a specific projection
    #fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=ccrs.PlateCarree())

    # Plot the data on a latitude and longitude scale
    im = ax.contourf(lon, lat, var, vmin=0, vmax=1, transform=ccrs.PlateCarree(), cmap='PuRd', levels=15)#, norm=matplotlib.colors.LogNorm(), levels=np.logspace(np.log10(var.min()), np.log10(var.max()), 10), extend='max')

    # Set the extent of the map to match your data
    ax.set_extent([-160, -70, -60, 0], crs=ccrs.PlateCarree())

    # Add parallels and meridians
    ax.gridlines(draw_labels=False, linewidth=0.5, color='grey', alpha=0.5, linestyle='-')
    ax.set_xticks(np.arange(-160, -60, 20), crs=ccrs.PlateCarree())
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
    if layer == "surf":
        layer_name = "Epipelagic"
    if layer == "depth":
        layer_name = "Mesopelagic"
    ax.set_title(rf'Turnover')
    #ax.set_xlabel(r'Longitude')
    #ax.set_ylabel(r'Latitude')

    return im, ax, cbar

def getis(plankton, layer, fig, rows, cols, pos):
    """
    :param plankton: "phyto" or "zoo"
    :param layer: "surf" or "depth"
    :return: fig, ax, cbar
    """

    if plankton == "phyto":
        if layer == "surf":
            path = "Data/getis_phyto_surf.nc"
            var_name = 'GiPhytosurface'
        elif layer == "depth":
            path = "Data/getis_phyto_depth.nc"
            var_name = 'GiPhytoDepth'
    elif plankton == "zoo":
        if layer == "surf":
            path = "Data/getis_zoo_surf.nc"
            var_name = 'Gizoosurface'
        elif layer == "depth":
            path = "Data/getis_zoo_depth.nc"
            var_name = 'GiZooDepth1.tif'

    data = nc.Dataset(path)

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
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=ccrs.PlateCarree())

    # Norm the data so 1 is neutral
    #norm=colors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=3.5)
    boundaries = np.array([-4, -3, -2, -1, 1, 2, 3, 4])
    norm = matplotlib.colors.BoundaryNorm(boundaries=boundaries, ncolors=10, extend='both')
    print(norm([-3.4, -2.5, 2, 0.4, -0.8]))

    # Plot the data on a latitude and longitude scale
    im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap='coolwarm', levels=[-3, -2, -1, 1, 2, 3], extend='both')#, norm=matplotlib.colors.LogNorm(), levels=np.logspace(np.log10(var.min()), np.log10(var.max()), 10), extend='max')

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
    ax.set_title(rf'Getis Ord')
    #ax.set_xlabel(r'Longitude')
    #ax.set_ylabel(r'Latitude')

    return im, ax, cbar


def SDM_spatial(plankton, layer, fig, rows, cols, pos=1):
    """
    :param plankton: "phyto" or "zoo"
    :param layer: "surf" or "depth"
    :return: fig, ax, cbar
    """
    path = f"Data/spatial{plankton}{layer}.nc"

    data = nc.Dataset(path)
    var_name = os.path.basename(path).split(".")[0] # Check this
    if var_name != 'spatialphytosurf':
        var_name = 'spatial'

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
    #fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=ccrs.PlateCarree())

    # Plot the data on a latitude and longitude scale
    im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap='plasma', norm=LogNorm(), levels=np.logspace(-4, np.log10(var.max()), 11), extend='min')

    # Set the extent of the map to match your data
    ax.set_extent([-160, -70, -60, 0], crs=ccrs.PlateCarree())

    # Add parallels and meridians
    ax.gridlines(draw_labels=False, linewidth=0.5, color='grey', alpha=0.5, linestyle='-')
    ax.set_xticks(np.arange(-160, -60, 20), crs=ccrs.PlateCarree())
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
    cbar.ax.yaxis.set_major_locator(plt.LogLocator())

    # Set the title and labels
    if layer == "surf":
        layer_name = "Epipelagic"
    if layer == "depth":
        layer_name = "Mesopelagic"
    ax.set_title(rf'Richness')
    #ax.set_xlabel(r'Longitude')
    #ax.set_ylabel(r'Latitude')

    return im, ax, cbar

def SDM_prediction(plankton, layer, fig, rows, cols, pos=1):
    """
    :param plankton: "phyto" or "zoo"
    :param layer: "surf" or "depth"
    :return: fig, ax, cbar
    """
    path = f"Data/prediction{plankton}{layer}.nc"

    data = nc.Dataset(path)
    var_name = os.path.basename(path).split(".")[0] # Check this
    if plankton =='zoo':
        var_name = 'prediction'
    elif var_name == 'predictionphytodepth':
        var_name = 'predictionphyto'

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
    #fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=ccrs.PlateCarree())

    # Plot the data on a latitude and longitude scale
    im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap='magma')

    # Set the extent of the map to match your data
    ax.set_extent([-160, -70, -60, 0], crs=ccrs.PlateCarree())

    # Add parallels and meridians
    ax.gridlines(draw_labels=False, linewidth=0.5, color='grey', alpha=0.5, linestyle='-')
    ax.set_xticks(np.arange(-160, -60, 20), crs=ccrs.PlateCarree())
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
    #cbar.ax.yaxis.set_major_locator(plt.LogLocator())

    # Set the title and labels
    if layer == "surf":
        layer_name = "Epipelagic"
    if layer == "depth":
        layer_name = "Mesopelagic"
    ax.set_title(rf'Richness')
    #ax.set_xlabel(r'Longitude')
    #ax.set_ylabel(r'Latitude')

    return im, ax, cbar
    """
    :param plankton: "phyto" or "zoo"
    :param layer: "surf" or "depth"
    :return: fig, ax, cbar
    """
    path = f"Data/prediction{plankton}{layer}.nc"

    data = nc.Dataset(path)
    var_name = os.path.basename(path).split(".")[0] # Check this
    if var_name != 'predictionphytosurf':
        var_name = 'prediction'

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
    #fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=ccrs.PlateCarree())

    # Plot the data on a latitude and longitude scale
    im = ax.contourf(lon, lat, var, transform=ccrs.PlateCarree(), cmap='magma')

    # Set the extent of the map to match your data
    ax.set_extent([-160, -70, -60, 0], crs=ccrs.PlateCarree())

    # Add parallels and meridians
    ax.gridlines(draw_labels=False, linewidth=0.5, color='grey', alpha=0.5, linestyle='-')
    ax.set_xticks(np.arange(-160, -60, 20), crs=ccrs.PlateCarree())
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
    cbar.ax.yaxis.set_major_locator(plt.LogLocator())

    # Set the title and labels
    if layer == "surf":
        layer_name = "Epipelagic"
    if layer == "depth":
        layer_name = "Mesopelagic"
    ax.set_title(rf'Richness')
    #ax.set_xlabel(r'Longitude')
    #ax.set_ylabel(r'Latitude')

    return im, ax, cbar
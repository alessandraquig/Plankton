import numpy as np
from netCDF4 import Dataset
import os

# Open the netCDF file
file_name = 'SST_1993-2022.nc'
nc = Dataset(f'Data/{file_name}', 'r')
var_name = os.path.basename(file_name).split("_1993-2022.nc")[0]
print(f"Variable name: {var_name}")


# Read the latitude, longitude, and time variables
lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:]
time = nc.variables['time'][:]

# Read the data variable
data = nc.variables['thetao_mean'][:]

# Calculate the temporal mean at each grid cell
mean_data = np.nanmean(data, axis=0)

# Check that it's giving the right thing
print(mean_data.shape)
# Create a new netCDF file
new_nc = Dataset(f'Data/{var_name}_mean.nc', 'w', format='NETCDF4')

# Create dimensions for latitude and longitude
new_nc.createDimension('latitude', len(lat))
new_nc.createDimension('longitude', len(lon))

# Create variables for latitude, longitude, and mean_data
lat_var = new_nc.createVariable('latitude', 'f4', ('latitude',))
lon_var = new_nc.createVariable('longitude', 'f4', ('longitude',))
data_var = new_nc.createVariable('temperature', 'f4', ('latitude', 'longitude'))

# Assign values to latitude, longitude, and mean_data variables
lat_var[:] = lat
lon_var[:] = lon
data_var[:] = mean_data

print(new_nc.shape)

# Close the new netCDF file
new_nc.close()

# Close the netCDF file
nc.close()
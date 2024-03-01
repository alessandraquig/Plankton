import netCDF4 as nc

# Open the NetCDF file
file_path = "Data/chlorophyll_mean.nc"
dataset = nc.Dataset(file_path)

# View the metadata
print("Metadata:")
print(dataset)

# View the first three lines of data
print("Data:")
for var_name in dataset.variables:
    var = dataset.variables[var_name]
    print(f"{var_name}:")
    print(var[:3])

# Close the NetCDF file
dataset.close()
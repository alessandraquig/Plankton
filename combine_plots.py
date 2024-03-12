import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import numpy as np
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import plotting_functions as pf
matplotlib.use('Agg')

from matplotlib import rc
this_rc_params = {
    "text.usetex": True,
    "font.family": "roman"
}
plt.rcParams.update(this_rc_params)

if __name__ == '__main__':
    
    # Plot environmental data (Figure 1)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    im_DO, axes[0,0], cbar_DO = pf.DO(layer='depth')
    im_SST, axes[0,1], cbar_SST = pf.SST()
    im_SiNO3, axes[1,0], cbar_SiNO3 = pf.SiNO3(layer='surf')
    im_chlorophyll, axes[1,1], cbar_chlorophyll = pf.chlorophyll()
    plt.tight_layout()
    plt.savefig('environmental_factors.pdf', type='pdf')
    plt.show()
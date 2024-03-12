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
    rows = 2
    cols = 2
    fig_enviro = plt.figure(figsize=(10, 8))
    fig_enviro.suptitle("Environmental Factors", fontsize=16)
    im_DO, ax_DO, cbar_DO = pf.DO(layer='depth', fig=fig_enviro, rows=rows, cols=cols, pos=1)
    ax_DO.set_title(r'Dissolved Oxygen')
    im_SST, ax_SST, cbar_SST = pf.SST(fig=fig_enviro, rows=rows, cols=cols, pos=2)
    ax_SST.set_title(r'Sea Surface Temperature ($^{\circ}$C)')
    im_SiNO3, ax_SiNO3, cbar_SiNO3 = pf.SiNO3(layer='surf', fig=fig_enviro, rows=rows, cols=cols, pos=3)
    ax_SiNO3.set_title(r'Silicate/Nitrate Ratio')
    im_chlorophyll, ax_chlorophyll, cbar_chlorophyll = pf.chlorophyll(fig=fig_enviro, rows=rows, cols=cols, pos=4)
    ax_chlorophyll.set_title(r'Chlorophyll Concentration ($mg/m^{3}$)')
    plt.tight_layout()
    plt.savefig('Output/environmental_factors.pdf')

    # Plot biological data (Figure 2)
    


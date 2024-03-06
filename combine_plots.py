import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import numpy as np
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import plot_turnover
import plot_bathymetry
import plot_chlorophyll
#import plot_DO
#import plot_getis
import plot_nestedness
import plot_richness
#import plot_SiNO3
import plot_SST
matplotlib.use('Agg')

from matplotlib import rc
this_rc_params = {
    "text.usetex": True,
    "font.family": "roman"
}
plt.rcParams.update(this_rc_params)

if __name__ == '__main__':
    layer = 'surf'
    plankton = 'phyto'

    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    im_rich, axes[0], cbar_rich = plot_richness.plot_richness(layer=layer, plankton=plankton)
    im_nest, axes[1], cbar_nest = plot_nestedness.plot_nestedness(layer=layer, plankton=plankton)
    plt.tight_layout()
    plt.savefig('output.png')
    plt.show()
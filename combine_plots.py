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
    ax_DO.set_title(r'Dissolved Oxygen ($mmol/m^{3}$)')
    im_SST, ax_SST, cbar_SST = pf.SST(fig=fig_enviro, rows=rows, cols=cols, pos=2)
    ax_SST.set_title(r'Sea Surface Temperature ($^{\circ}$C)')
    im_SiNO3, ax_SiNO3, cbar_SiNO3 = pf.SiNO3(layer='surf', fig=fig_enviro, rows=rows, cols=cols, pos=3)
    ax_SiNO3.set_title(r'Silicate/Nitrate Ratio')
    im_chlorophyll, ax_chlorophyll, cbar_chlorophyll = pf.chlorophyll(fig=fig_enviro, rows=rows, cols=cols, pos=4)
    ax_chlorophyll.set_title(r'Chlorophyll Concentration ($mg/m^{3}$)')
    plt.tight_layout()
    plt.savefig('Output/environmental_factors.pdf')

    # Plot biological data (Figure 2)
    rows = 2
    cols = 2
    layer = 'surf'
    plankton = 'phyto'
    fig_bio = plt.figure(figsize=(10, 8))
    fig_bio.suptitle(r"Epipelagic Phytoplankton", fontsize=16)
    im_rich, ax_rich, cbar_rich = pf.richness(plankton=plankton, layer=layer, fig=fig_enviro, rows=rows, cols=cols, pos=1)
    ax_rich.set_title(r'Richness')
    im_getis, ax_getis, cbar_getis = pf.getis(plankton=plankton, layer=layer, fig=fig_enviro, rows=rows, cols=cols, pos=2)
    ax_getis.set_title(r'Gi\*')
    im_nest, ax_nest, cbar_nest = pf.nestedness(plankton=plankton, layer=layer, fig=fig_enviro, rows=rows, cols=cols, pos=3)
    ax_nest.set_title(r'Nestedness')
    im_turn, ax_turn, cbar_turn = pf.turnover(plankton=plankton, layer=layer, fig=fig_enviro, rows=rows, cols=cols, pos=4)
    ax_turn.set_title(r'Turnover')
    plt.tight_layout()
    plt.savefig('Output/epipelagic_phytoplankton.pdf')

    rows = 1
    cols = 1
    layer = 'surf'
    plankton = 'phyto'
    fig_bio = plt.figure(figsize=(10, 8))
    fig_bio.suptitle(r"Epipelagic Phytoplankton", fontsize=16)
    im_getis, ax_getis, cbar_getis = pf.getis(plankton=plankton, layer=layer, fig=fig_enviro, rows=rows, cols=cols, pos=1)
    ax_getis.set_title(r'Gi\*')
    plt.savefig('Output/epipelagic_phytoplankton.pdf')


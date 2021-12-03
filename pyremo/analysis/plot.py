

def plot_seasons(da, vmin=None, vmax=None, extent=None, cmap='bwr', transform=None,
                projection=None, crs_extent=None, borders=False, xlocs=range(-180,180,10), ylocs=range(-90,90,10),
                aspect='auto'):
    """plot seasonal means"""
    from matplotlib import pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cf
    if projection is None:
        projection=ccrs.PlateCarree()
    if transform is None:
        transform=ccrs.PlateCarree()
    if crs_extent is None:
        crs_extent = transform
    #plt.subplots_adjust(hspace=1.5, wspace=1.0)
    fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(
        ncols=2, nrows=2, subplot_kw={'projection': projection}, 
        figsize=(18,14))
    
    axes = (ax1, ax2, ax3, ax4)
    for ax in axes:
        #ax.set_axis_off()
        if extent is None:
            ax.set_extent([da.rlon.min(), da.rlon.max(), da.rlat.min(), da.rlat.max()], crs=transform)
        else:
            ax.set_extent(extent)
        ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', 
                      xlocs=xlocs, ylocs=ylocs)
        ax.coastlines(resolution='50m', color='black', linewidth=1)
        if aspect is not None:
            ax.set_aspect(aspect)
        if borders: ax.add_feature(cf.BORDERS)
    for season, ax in zip(da.season, axes):
        im = da.sel(season=season).plot(ax=ax, vmin=vmin, vmax=vmax, cmap=cmap, 
                                        transform=transform, add_colorbar=False)
    cbar = fig.colorbar(im, ax=axes)
    return plt
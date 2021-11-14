



def plot_seasons(da, pole, vmin=None, vmax=None, cmap='coolwarm'):
    """plot seasonal means"""
    from matplotlib import pyplot as plt
    import cartopy.crs as ccrs
    transform = ccrs.RotatedPole(pole[0], pole[1])
    projection = transform
    #plt.subplots_adjust(hspace=1.5, wspace=1.0)
    fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(
        ncols=2, nrows=2, subplot_kw={'projection': projection}, 
        figsize=(18,14))
    
    axes = (ax1, ax2, ax3, ax4)
    for ax in axes:
        #ax.set_axis_off()
        #ax.set_extent([da.rlon.min(), da.rlon.max(), da.rlat.min(), da.rlat.max()], crs=transform)
        ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', 
                      xlocs=range(-180,180,10), ylocs=range(-90,90,10))
        ax.coastlines(resolution='110m', color='black', linewidth=1)
    for season, ax in zip(da.season, axes):
        im = da.sel(season=season).plot(ax=ax, vmin=vmin, vmax=vmax, cmap=cmap, 
                                        transform=transform, add_colorbar=False)
    cbar = fig.colorbar(im, ax=axes)
    return plt
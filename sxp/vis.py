import numpy as np
import matplotlib.pyplot as pl
from matplotlib import animation as ani
from matplotlib.colors import LogNorm

from photutils.centroids import centroid_com
from photutils.centroids import centroid_1dg
from photutils.centroids import centroid_2dg

from tqdm import tqdm


def movie(frames, name, fps=60, dpi=100):

    FFMpegWriter = ani.writers['mencoder']
    writer = FFMpegWriter(fps=fps)
    fig = pl.figure(frameon=False)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    im = ax.imshow(frames[0], interpolation='none',
        cmap=pl.cm.gray, norm=LogNorm())

    with writer.saving(fig, name, dpi):
        for frame in tqdm(frames):
            try:
                im.set_data(frame)
                writer.grab_frame()
            except:
                pass

def plot_stacked(cube, outpath=None, dpi=240):

    med_im = np.median(cube, axis=0)

    def fix_nan(im):
        im[np.isnan(im)] = np.median(im[~np.isnan(im)])
    fix_nan(med_im)

    pl.figure(figsize=(10,10))
    pl.imshow(med_im, cmap=pl.cm.gray, interpolation='none', norm=LogNorm())
    xl, yl = pl.xlim(), pl.ylim()

    mask = np.ones_like(med_im).astype(bool)
    x0 = y0 = 10
    x1 = y1 = 21
    mask[x0:x1,y0:y1] = False

    cens = []
    ma = np.ma.masked_array(med_im, mask)
    masked = ma.filled(0)
    s = 20
    pl.plot(15, 15, 'yo', fillstyle='none', ms=30, markeredgewidth=3)
    for centroid in centroid_1dg, centroid_2dg, centroid_com:
        # cen = centroid(med_im, mask=mask)
        cen = centroid(masked)
        cens.append(cen)
        x, y = cen
        print x, y
        pl.plot(x, y, 'o', ms=s, fillstyle='none',
            markeredgewidth=3, label=str(centroid).split()[1])
        s -= 5
    cen_med = np.nanmedian(cens, axis=0).tolist()
    pl.legend(numpoints=1)
    pl.xlim(*xl), pl.ylim(*yl)
    pl.tight_layout()
    if outpath:
        pl.savefig(outpath, dpi=dpi)
        pl.close()
    else:
        pl.show()
    return cen_med


def centroids(t, x, y, outpath=None, ms=2):

    t = np.linspace(0, (t[-1] - t[0]) * 24, t.size)

    pl.subplot(211)
    pl.plot(t, x, 'k.', ms=ms)
    pl.ylabel('x centroid')
    pl.gca().xaxis.set_ticklabels(['']*len(pl.gca().xaxis.get_ticklabels()))
    pl.subplot(212)
    pl.plot(t, y, 'k.', ms=ms)
    pl.xlabel('time')
    pl.ylabel('y centroid')
    pl.tight_layout()

    if outpath:
        pl.savefig(outpath)
        pl.close()
    else:
        pl.show()


def multi_ts(t, ts_list, outpath=None, ms=2):

    colors = [pl.cm.gist_earth(i) for i in np.linspace(0,1,len(ts_list))]

    for i,f in enumerate(ts_list):
        pl.plot(t, f, '.', color=colors[i], ms=ms)
    pl.xlabel('time')
    pl.ylabel('flux')

    if outpath:
        pl.savefig(outpath)
        pl.close()
    else:
        pl.show()


def multi_ts2(time, fluxes_r, outpath=None):

    nr = len(fluxes_r)

    colors = [pl.cm.RdBu_r(i) for i in np.linspace(0, 1, nr)]
    alphas = np.linspace(0.25, 0.75, nr)

    sizes = range(2, nr+2)
    sizes.reverse()

    pl.figure(figsize=(15,10))
    for i,flux in enumerate(fluxes_r):
        pl.plot(time, flux, 'o', color=colors[i], ms=sizes[i], alpha=alphas[i])
    pl.xlabel('time')
    pl.ylabel('flux')

    if outpath:
        pl.savefig(outpath)
        pl.close()
    else:
        pl.show()



import os
import sys
import shutil
import pickle
import numpy as np
import pandas as pd
from photutils.centroids import centroid_2dg
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as pl



def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if shutil.fnmatch.fnmatch(basename, pattern):
                filepath = os.path.join(root, basename)
                yield filepath


def ignore_oserror(f):
    """
    rather silly decorator for ignoring OSError exceptions
    """
    def wrapper(*args):
        try:
            f(*args)
        except OSError as e:
            print(e)
    return wrapper


@ignore_oserror
def mkdir(d):
    os.makedirs(d)


def mkdirs(*args):
    for d in args:
        mkdir(d)


def load_data(data_dir, transit_aor):

    """
    for reading the current output format of sxp photometry.
    """

    d = pickle.load(open('{}/{}_phot.pkl'.format(data_dir, transit_aor), 'rb'))
    return [d[i] for i in 'cube, time, flux, radii, unc, cen'.split(', ')]


def get_pix(cube, cen=None, geom='3x3', normalize=True):

    """
    assumes the cube has shape (n,k,k), where n is the number of frames and k
    is the width of each (square) frame, and that k > 5.
    """

    def fix_nan(im):
        im[np.isnan(im)] = np.median(im[~np.isnan(im)])

    if not cen:
        med_im = np.median(cube, axis=0)
        fix_nan(med_im)
        cx, cy = centroid_2dg(med_im)
        cx, cy = list(map(int, list(map(round, [cx, cy]))))
        print("centroid: {}, {}".format(cx, cy))
    else:
        cx, cy = list(map(int, cen))
    if geom == '3x3':
        i = 1
    elif geom == '5x5':
        i = 2
    else:
        raise ValueError("geometry not supported")

    x0 = cx - i
    x1 = cx + i + 1
    y0 = cy - i
    y1 = cy + i + 1
    sub_cube = cube[:, x0:x1, y0:y1]
    pixels = sub_cube.reshape(cube.shape[0], sub_cube.shape[1]**2)

    if normalize:
        for i in pixels: i /= i.sum()

    return pixels


def df_from_pickle(picklefile, radius, pix=False, geom='3x3', normalize=True):

    d = pickle.load(open(picklefile, 'rb'))
    idx = d['radii'].index(radius)
    t = d['time']
    f = d['flux'][idx]
    s = d['unc'][idx]
    df = pd.DataFrame(dict(t=t,f=f,s=s))
    df = df[['t', 'f', 's']]
    if not pix:
        return df
    cube = d['cube']
    pix = get_pix(cube, geom=geom, normalize=normalize)
    for i in range(pix.shape[1]):
        key = 'p{}'.format(i)
        df[key] = pix[:,i]
    return df


def to_mags(flux, zp):

    """
    Converts flux to magnitudes using the input zero point,
    which must be in the same units as flux (i.e. Janskies)
    """

    return -2.5 * np.log10(flux/zp)


def spz_jy_to_mags(jy, ch):

    """
    Converts IRAC ch1 and ch2 flux in Janskies to magnitudes.
    """

    if ch==1:
        zp = 280.9
    elif ch==2:
        zp = 179.7
    else:
        raise ValueError('ch must be either 1 or 2')
    return to_mags(jy,zp)


def binned(a, binsize, fun=np.mean):
    return np.array([fun(a[i:i+binsize], axis=0) \
        for i in range(0, a.shape[0], binsize)])


def beta(flux, timestep, start_min=5, stop_min=20):

    """
    flux : measured flux (presumably normalized)
    timestep : time interval between datapoints in seconds
    """

    assert timestep < start_min * 60
    ndata = len(flux)

    sigma1 = np.std(flux)

    min_bs = int(start_min * 60 / timestep)
    max_bs = int(stop_min * 60 / timestep)

    betas = []
    for bs in range(min_bs, max_bs + 1):
        nbins = int(round(ndata / bs))
        if nbins < 2:
            continue
        sigmaN_theory = sigma1 / np.sqrt(bs) * np.sqrt( nbins / (nbins - 1) )
        sigmaN_actual = np.std(binned(flux, bs))
        beta = sigmaN_actual / sigmaN_theory
        betas.append(beta)

    return np.median(betas)


def compute_best_radius(pkl, n_seg=10, fp=None):

    if fp is not None: fig, axs = pl.subplots(1, 2, figsize=(10,3))

    t = pkl[b'time']
    flux = pkl[b'flux']
    rad = pkl[b'radii']
    timestep = np.diff(t).mean() * 86400

    r_std, r_beta =[], []

    for i in range(n_seg):
        f_std = []
        f_beta = []

        n = int(len(t)/n_seg)
        for r in rad:
            idx = rad.index(r)
            f_std.append(flux[idx][i*n:(i+1)*n].std())
            f_beta.append(beta(flux[idx][i*n:(i+1)*n], timestep))

        r_std_opt = rad[f_std.index(min(f_std))]
        f_std_opt = f_std[f_std.index(min(f_std))]
        if fp is not None:
            axs[0].plot(rad, f_std)
            axs[0].plot(r_std_opt, f_std_opt, 'ko')
            axs[0].set_title(r'$\sigma$')

        r_beta_opt = rad[f_beta.index(min(f_beta))]
        f_beta_opt = f_beta[f_beta.index(min(f_beta))]
        if fp is not None:
            axs[1].plot(rad, f_beta)
            axs[1].plot(r_beta_opt, f_beta_opt, 'ko')
            axs[1].set_title(r'$\beta$')

        r_std.append(r_std_opt)
        r_beta.append(r_beta_opt)

    if fp is not None:
        axs[0].vlines(np.median(r_std), *axs[0].get_ylim(), linestyle='--', color='gray')
        axs[1].vlines(np.median(r_beta), *axs[1].get_ylim(), linestyle='--', color='gray')

    if fp is not None:
        pl.savefig(fp)
        pl.close()

    r_opt_std = round(np.median(r_std),1)
    r_opt_beta = round(np.median(r_beta),1)

    print(("optimal for white noise: {}".format(r_opt_std)))
    print(("optimal for red noise: {}".format(r_opt_beta)))

    idx = rad.index(round(np.mean([r_opt_std,r_opt_beta]),1))
    print(("using radius: {}".format(rad[idx])))
    f = flux[idx]

    return t, f

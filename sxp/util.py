import os
import sys
import shutil
import pickle
import numpy as np
from photutils.morphology import centroid_com


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


def get_pix(cube, cen=None, geom='3x3', normalize=False):

    """
    assumes the cube has shape (n,k,k), where n is the number of frames and k
    is the width of each (square) frame, and that k > 5.
    """

    def fix_nan(im):
        im[np.isnan(im)] = np.median(im[~np.isnan(im)])

    if not cen:
        med_im = np.median(cube, axis=0)
        fix_nan(med_im)
        cx, cy = centroid_com(med_im)
        cx, cy = map(int, map(round, [cx, cy]))
        print "centroid: {}, {}".format(cx, cy)
    else:
        cx, cy = map(int, cen)
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

import numpy as np
from scipy import stats
from astropy.io import fits
# from astropy import units as u
from astropy.stats import sigma_clip
from photutils import aperture_photometry, CircularAperture
from photutils.morphology import centroid_com, centroid_1dg, centroid_2dg
from tqdm import tqdm
from util import find_files


def read_data(data_dir):

	search_str = "*I2*_bcd.fits"
	files = list(find_files(data_dir, search_str))

	hdr = fits.getheader(files[0]).copy()
	subarray = True
	if hdr['NAXIS1'] == 256:
		subarray = False

	# read data and header keyword values
	print "reading data"
	data, time_bmjd, obs_length_s = [], [], []
	for fp in tqdm(files):
		hdr = fits.getheader(fp).copy()
		im = fits.getdata(fp).copy()
		data.append(im)
		time_bmjd.append(hdr['BMJD_OBS'])
		obs_length_s.append(hdr['ATIMEEND'] - hdr['AINTBEG'])

	# make into arrays and sort by start time
	data, time_bmjd, obs_length_s = map(np.array, [data, time_bmjd, obs_length_s])
	time_bjd = time_bmjd + 2400000.5
	idx = np.argsort(time_bmjd)
	time_bmjd = time_bmjd[idx]
	time_bjd = time_bjd[idx]
	obs_length_s = obs_length_s[idx]
	data = data[idx]
	assert np.all(np.argsort(time_bmjd) == np.arange(idx.size))

	# calculate times for each individual frame
	time = []
	for tup in zip(time_bjd, obs_length_s):
		t, length_s = tup
		length_d = length_s * (1./60) * (1./60) * (1./24)
		if subarray:
			frame_length = length_d / 64.
			for i in range(64):
				time.append(t + frame_length * i)
		else:
			time.append(t + length_d / 2.)

	time = np.array(time)
	assert all([time[i] < time[i+1] for i in range(time.size-1)]) # another paranoid check

	# create data cube
	if subarray:
		cube = data.reshape(data.shape[0] * data.shape[1], data.shape[2], data.shape[3])
		pix_outliers(cube)
	else:
		cube = np.array([i[-40:-8,8:40] for i in data])
	return time, cube


def pix_outliers(cube, step=64, sigma_level=3, final_iter=True, verbose=False):

    nframes = cube.shape[0]
    if verbose:
		print "num frames: {}".format(nframes)

    for i in tqdm(range(nframes/step)):

        if (i+1)*step <= nframes:
            chunk = cube[i*step:(i+1)*step]
            if verbose:
				print "processing frames {} - {}".format(i*step, (i+1)*step-1)
        else:
            chunk = cube[i*step:]
            if verbose:
				print "processing frames {} - {}".format(i*step, nframes-1)

        sd = np.std(chunk, axis=0)
        md = np.median(chunk, axis=0)

        for frame in chunk:
            idx = np.abs(frame - md) > sigma_level * sd
            frame[idx] = md[idx]

	if final_iter:
		cubec = sigma_clip(cube, sigma=5, axis=0)
		mask = cubec.mask

		med_im = np.nanmedian(cube, axis=0)

		for i,m in enumerate(mask):
		    cube[i][m] = med_im[m]


def get_centroids(cube, centroid=centroid_com):

	im = cube[0]
	x0 = y0 = 10
	x1 = y1 = 21

	mask = np.ones_like(im).astype(bool)
	mask[x0:x1,y0:y1] = False

	centroids = []
	for im in tqdm(cube):
		ma = np.ma.masked_array(im, mask)
		masked = ma.filled(0)
		centroids.append(centroid(masked))

	return np.array(centroids)


def centroid_sigclip(time, cube, centroid, final_iter=True, verbose=True):

	# mask bad centroids, if any exist
	idx = np.array([np.isnan(i).any() for i in centroid])
	time = time[~idx]
	cube = cube[~idx]
	centroid = centroid[~idx]

	# number of minutes per time step
	time_step_m = np.diff(time).mean() * 24 * 60

	# window size: number of frames per 40 min
	# (40 min is the period of the battery heater effect)
	ws = int(40 / time_step_m)

	nframes = cube.shape[0]
	nsteps = cube.shape[0] / ws

	assert time.size == centroid.shape[0]
	flagged = np.zeros_like(time).astype(bool)

	x, y = zip(*centroid)
	print "removing centroid outliers"
	for i in tqdm(range(nsteps)):
	    if (i+1) * ws <= nframes:
	        x_chunk = x[i*ws:(i+1)*ws]
	        y_chunk = x[i*ws:(i+1)*ws]
	        flag_chunk = flagged[i*ws:(i+1)*ws]
	    else:
	        x_chunk = y[i*ws:]
	        y_chunk = y[i*ws:]
	        flag_chunk = flagged[i*ws:]

	    x_median = np.median(x_chunk)
	    x_sigma = np.std(x_chunk)
	    y_median = np.median(y_chunk)
	    y_sigma = np.std(y_chunk)

	    idx = (np.abs(x_chunk - x_median) > 3 * x_sigma) | (np.abs(y_chunk - y_median) > 3 * y_sigma)
	    flag_chunk[idx] = True

	if verbose:
	    print "{} centroid outliers removed".format(flagged.sum())

	time = time[~flagged]
	cube = cube[~flagged]
	centroid = centroid[~flagged]

	if final_iter:
		x, y = zip(*centroid)
		xc = sigma_clip(x, sigma=5)
		yc = sigma_clip(y, sigma=5)
		mask = xc.mask | yc.mask
		time = time[~mask]
		cube = cube[~mask]
		centroid = centroid[~mask]
		if verbose:
		    print "{} centroid outliers removed".format(mask.sum())

	return time, cube, centroid


def multi_aperture(radii, cube, centroids):

	fluxes_r = []

	print "computing photometry"
	for r in tqdm(radii):

		flux = []
		for tup in zip(cube, centroids):
			im, cen = tup
			apertures = CircularAperture(cen, r=r)
			phot_table = aperture_photometry(im, apertures)
			converted_aperture_sum = phot_table['aperture_sum']
			flux.append(converted_aperture_sum)
		flux = np.array(flux)
		flux = flux.reshape(flux.shape[0])
		fluxes_r.append(flux)

	return fluxes_r


def flux_sigmaclip(time, cube, centroid, radii, fluxes_r, unc_r,
	final_iter=True, verbose=True):

	# number of minutes per time step
	time_step_m = np.diff(time).mean() * 24 * 60
	ws = int(40 / time_step_m)
	nframes = cube.shape[0]
	nsteps = cube.shape[0] / ws

	flux_masks = []
	for i in range(len(fluxes_r)):

	    flux = fluxes_r[i]
	    flagged = np.zeros_like(flux).astype(bool)

	    for j in range(nsteps):
	        if (j+1) * ws <= nframes:
	            chunk = flux[j*ws:(j+1)*ws]
	            flag_chunk = flagged[j*ws:(j+1)*ws]
	        else:
	            chunk = flux[j*ws:]
	            flag_chunk = flagged[j*ws:]

	        f_median = np.median(chunk)
	        f_sigma = np.std(chunk)

	        idx = np.abs(chunk - f_median) > 3 * f_sigma
	        flag_chunk[idx] = True

		flux_masks.append(flagged)

	if verbose:
		for i in range(len(fluxes_r)):
			msg = "{} flux outliers for radius {}"
			print msg.format(flux_masks[i].sum(), radii[i])

	# TODO: currently throwing out a time step if *ANY* of the
	# apertures produced an outlying flux measurement
	flagged = np.array(flux_masks).any(axis=0)
	cube = cube[~flagged]
	time = time[~flagged]
	for i in range(len(fluxes_r)):
		fluxes_r[i] = fluxes_r[i][~flagged]
		unc_r[i] = unc_r[i][~flagged]
	centroid = centroid[~flagged]

	if final_iter:
		masks = []
		for flux in fluxes_r:
			fluxc = sigma_clip(flux, sigma=5)
			masks.append(fluxc.mask)
		flagged = np.array(masks).any(axis=0)
		cube = cube[~flagged]
		time = time[~flagged]
		for i in range(len(fluxes_r)):
			fluxes_r[i] = fluxes_r[i][~flagged]
			unc_r[i] = unc_r[i][~flagged]
		centroid = centroid[~flagged]
		if verbose:
			for i in range(len(fluxes_r)):
				msg = "{} flux outliers for radius {}"
				print msg.format(masks[i].sum(), radii[i])

	return time, cube, centroid, fluxes_r, unc_r


def subtract_bg(cube):

	mask = np.zeros(cube[0].shape).astype(bool)
	mask[6:-6, 6:-6] = True
	mask[12:16,:] = True
	mask[:,13:15] = True
	mask[:,31] = True
	unc = []
	print "subtracting background, computing uncertainties, normalizing flux"
	for im in tqdm(cube):
		good = ~np.isnan(im) & ~mask
		# knutson: iterative clipping of 3-sigma outliers prior to gaussian fit
		clipped, lo, up = stats.sigmaclip(im[good], low=3., high=3.)
		fit = stats.norm.fit(clipped)
		unc.append(fit[1])
		im -= fit[0]
	unc = np.array(unc)
	return cube, unc


def get_unc(unc, radii, fluxes_r):

	unc_r = []

	for i in range(len(radii)):

		r = radii[i]
		unc_i = []
		for j in range(len(fluxes_r[0])):
			# multiply uncertainty/pixel by number of pixels in aperture
			u = unc[j] * np.pi * r ** 2
			f = fluxes_r[i][j]
			unc_i.append(u/f)

		unc_r.append(np.array(unc_i))

		# now normalize the flux
		fluxes_r[i] /= np.median(fluxes_r[i])

	return fluxes_r, unc_r
from spectral_cube import SpectralCube
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
from scipy.optimize import curve_fit
from scipy import *
import time
import pprocess
from astropy.convolution import convolve
import radio_beam

def gauss_fitter(region = 'Cepheus_L1251', direct = '/Users/jkeown/Desktop/GAS_dendro/', SN_thresh = 3.0, mol = 'C2S', peak_channels = [222,270], convolve=False, use_old_conv=False, parallel = 3):

    	"""
    	Fit a Gaussian to non-NH3 emission lines from GAS.
    	It creates a cube for the best-fit Gaussian, a cube 
    	for the best-fit Gaussian with noise added back into 
    	the spectrum, and a parameter map of Tpeak, Vlsr, and FWHM
    
    	Parameters
    	----------
    	region : str
        	Name of region to reduce
    	direct : str
        	Directory where emission cubes are stored
		Also the location where new maps are saved
    	SN_thresh : float
        	Lowest signal-to-noise pixels to include in the line-fitting
    	mol : str
        	name of molecule to fit
    	peak_channels : [int, int]
        	Range of channels over which emission peaks may be found 
    	convolve : bool or float
        	If not False, specifies the beam-size to convolve the original map with
		Beam-size must be given in arcseconds
    	use_old_conv : bool
        	If True, use an already convolved map with name:
		region + '_' + mol + '_conv.fits'
		This convolved map must be in units of km/s
    	parallel : int
		Maximum number of simultaneous processes desired 
    	"""

	# Load the spectral cube and convert to velocity units
	cube = SpectralCube.read(direct + region + '_' + mol + '_base_DR2.fits')
	cube_km = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio')

	# If desired, convolve map with larger beam 
	# or load previously created convolved cube
	if convolve!=False:
		cube = SpectralCube.read(direct + region + '_' + mol + '_base_DR2.fits')
		cube_km_1 = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio')
		beam = radio_beam.Beam(major=convolve*u.arcsec, minor=convolve*u.arcsec, pa=0*u.deg)
		cube_km = cube_km_1.convolve_to(beam)
		cube_km.write(direct + region + '_' + mol + '_conv.fits', format='fits', overwrite=True)
	if use_old_conv!=False:
		cube_km = SpectralCube.read(direct + region + '_' + mol + '_conv.fits')
	
	# Define the spectral axis in km/s
	spectra_x_axis_kms = np.array(cube_km.spectral_axis) 

	# Create cubes for storing the fitted Gaussian profiles
	# and the Gaussians with noise added back into the spectrum
	header = cube_km.header
	cube_gauss = np.array(cube_km.unmasked_data[:,:,:])
	cube_gauss_noise = np.array(cube_km.unmasked_data[:,:,:])
	shape = np.shape(cube_gauss)

	# Set up a cube for storing fitted parameters
	param_cube = cube_gauss[0]
	param_cube = param_cube.reshape((1,) + param_cube.shape)
	param_cube = np.concatenate((param_cube, param_cube, param_cube, param_cube, param_cube, 	param_cube), axis=0)
	param_header = cube_km.header

	# Define the Gaussian profile
	def p_eval(x, a, x0, sigma):
    		return a*np.exp(-(x-x0)**2/(2*sigma**2))

	# Create some arrays full of NANs
	# To be used in output cubes if fits fail
	nan_array=np.empty(shape[0]) # For gauss cubes
	nan_array[:] = np.NAN
	nan_array2=np.empty(param_cube.shape[0]) # For param cubes
	nan_array2[:] = np.NAN

	# Loop through each pixel and find those
	# with SNR above SN_thresh
	x = []
	y = []
	pixels = 0
	for (i,j), value in np.ndenumerate(cube_gauss[0]):
     		spectra=np.array(cube_km.unmasked_data[:,i,j])
     		rms = np.std(np.append(spectra[0:(peak_channels[0]-1)], spectra[(peak_channels[1]+1):len(spectra)]))
     		if (max(spectra[peak_channels[0]:peak_channels[1]]) / rms) > SN_thresh:
            		pixels+=1
	    		x.append(i)
	    		y.append(j)
     		else:	
	    		cube_gauss[:,i,j]=nan_array
	    		param_cube[:,i,j]=nan_array2
			cube_gauss_noise[:,i,j]=nan_array
	print str(pixels) + ' Pixels above SNR=' + str(SN_thresh) 

	# Define a Gaussian fitting function for each pixel
	# i, j are the x,y coordinates of the pixel being fit
	def pix_fit(i,j):
		spectra = np.array(cube_km.unmasked_data[:,i,j])
		# Use the peak brightness Temp within specified channel 
		# range as the initial guess for Gaussian height
		Tpeak = max(spectra[peak_channels[0]:peak_channels[1]])
		# Use the velocity of the brightness Temp peak as 
		# initial guess for Gaussian mean
		vpeak = spectra_x_axis_kms[peak_channels[0]:peak_channels[1]][np.where(spectra[peak_channels[0]:peak_channels[1]]==Tpeak)]
		rms = np.std(np.append(spectra[0:(peak_channels[0]-1)], spectra[(peak_channels[1]+1):len(spectra)]))
		err1 = np.zeros(shape[0])+rms
		# Create a noise spectrum based on rms of off-line channels
		# This will be added to best-fit Gaussian to obtain a noisy Gaussian 
		noise=np.random.uniform(-1.*rms,rms,len(spectra_x_axis_kms))
		# Define initial guesses for Gaussian fit
		guess = [Tpeak, vpeak, 0.3] # [height, mean, sigma]
		try:
			coeffs, covar_mat = curve_fit(p_eval, xdata=spectra_x_axis_kms, ydata=spectra, p0=guess, sigma=err1, maxfev=500)
			gauss = np.array(p_eval(spectra_x_axis_kms,coeffs[0], coeffs[1], coeffs[2]))
			noisy_gauss = np.array(p_eval(spectra_x_axis_kms,coeffs[0], coeffs[1], coeffs[2]))+noise
			params = np.append(coeffs, (covar_mat[0][0]**0.5, covar_mat[1][1]**0.5, covar_mat[2][2]**0.5)) 
			# params = ['Tpeak', 'VLSR','sigma','Tpeak_err','VLSR_err','sigma_err']

			# Don't accept fit if fitted parameters are non-physical or too uncertain
			if (params[0] < 0.01) or (params[3] > 1.0) or (params[2] < 0.05) or (params[5] > 0.5) or (params[4] > 0.75):
				noisy_gauss = nan_array
				gauss = nan_array
				params = nan_array2

			# Don't accept fit if the SNR for fitted spectrum is less than SNR threshold
			#if max(gauss)/rms < SN_thresh:
			#	noisy_gauss = nan_array
			#	gauss = nan_array
			#	params = nan_array2
		
		except RuntimeError:
			noisy_gauss = nan_array
			gauss = nan_array
			params = nan_array2
	
		return i, j, gauss, params, noisy_gauss

	# Parallel computation:
	nproc = parallel  # maximum number of simultaneous processes desired
	queue = pprocess.Queue(limit=nproc)
	calc = queue.manage(pprocess.MakeParallel(pix_fit))
	tic=time.time()
	counter = 0

	# Uncomment to see some plots of the fitted spectra
	#for i,j in zip(x,y):
		#pix_fit(i,j)
		#plt.plot(spectra_x_axis_kms, spectra, color='blue', drawstyle='steps')
		#plt.plot(spectra_x_axis_kms, gauss, color='red')
		#plt.show()
		#plt.close()

	# Begin parallel computations
	# Store the best-fit Gaussians and parameters 
	# in their correct positions in the previously created cubes
	for i,j in zip(x,y):
		calc(i,j)
	for i,j,gauss_spec,parameters,noisy_gauss_spec in queue:
		cube_gauss[:,i,j]=gauss_spec
		param_cube[:,i,j]=parameters
		cube_gauss_noise[:,i,j]=noisy_gauss_spec
		counter+=1
		print str(counter) + ' of ' + str(pixels) + ' pixels completed'
	print "%f s for parallel computation." % (time.time() - tic)

	# Save final cubes
	cube_final_gauss = SpectralCube(data=cube_gauss, wcs=cube_km.wcs, header=cube_km.header)
	cube_final_gauss.write(direct + 'gauss_cube_' + mol + '.fits', format='fits', overwrite=True)
	cube_final_gauss_noise = SpectralCube(data=cube_gauss_noise, wcs=cube_km.wcs, header=cube_km.header)
	cube_final_gauss_noise.write(direct + 'gauss_cube_noise_' + mol + '.fits', format='fits', overwrite=True)

	# Construct appropriate header for param_cube 
	param_header['NAXIS3'] = len(nan_array2)
	param_header['WCSAXES'] = 3
	param_header['CRPIX3'] = 1
	param_header['CDELT3'] = 1
	param_header['CRVAL3'] = 0
	param_header['PLANE1'] = 'Tpeak'
	param_header['PLANE2'] = 'VLSR'
	param_header['PLANE3'] = 'sigma'
	param_header['PLANE5'] = 'Tpeak_err'
	param_header['PLANE6'] = 'VLSR_err'
	param_header['PLANE7'] = 'sigma_err'

	cube_final_param = SpectralCube(data=param_cube, wcs=cube_km.wcs, header=param_header)
	cube_final_param.write(direct + 'param_cube_' + mol + '.fits', format='fits', overwrite=True)
	fits.writeto(direct + 'param_cube_' + mol + '.fits', param_cube, header=param_header, clobber=True)

### Examples ###
# Fit the HC5N data in Cepheus_L1251, without convolution
#gauss_fitter(region = 'Cepheus_L1251', direct = '/Users/jkeown/Desktop/GAS_dendro/', SN_thresh = 3.0, mol = 'HC5N', peak_channels = [402,460], convolve=False, use_old_conv=False)

# Convolve the HC5N data in Cepheus_L1251 to a spatial resolution of 64 arcseconds,
# then fit a Gaussian to all pixels above SNR=3
#gauss_fitter(region = 'Cepheus_L1251', direct = '/Users/jkeown/Desktop/GAS_dendro/', SN_thresh = 3.0, mol = 'HC5N', peak_channels = [402,460], convolve=64., use_old_conv=False)


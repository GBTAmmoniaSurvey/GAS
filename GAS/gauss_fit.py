from spectral_cube import SpectralCube
import pyspeckit
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
import sys
import os
import glob


def run_gauss_fits_all(file_extension='all_rebase3'):
        """
        Run the Gaussian fits on all observed regions
        """
        mol_list = ['C2S','HC5N','HC7N_22_21','HC7N_21_20','NH3_33']
        for mol in mol_list:
                print mol
                gauss_fitter(region="B1",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="B1E",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="B18",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="B59",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="Cepheus_L1228",mol=mol,vmin=-12,vmax=0,file_extension=file_extension)
                gauss_fitter(region="Cepheus_L1251",mol=mol,vmin=-12,vmax=0,file_extension=file_extension)
                gauss_fitter(region="CrAeast",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="CrAwest",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="HC2",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="IC348",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="IC5146",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="L1448",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="L1451",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="L1455",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="L1688",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="L1689",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="L1712",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="NGC1333",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="OrionA",mol=mol,vmin=0,vmax=15,file_extension=file_extension)
                gauss_fitter(region="OrionA_S",mol=mol,vmin=0,vmax=15,file_extension=file_extension)
                gauss_fitter(region="OrionB_NGC2023-2024",mol=mol,vmin=0,vmax=15,file_extension=file_extension)
                gauss_fitter(region="OrionB_NGC2068-2071",mol=mol,vmin=0,vmax=15,file_extension=file_extension)
                gauss_fitter(region="Perseus",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="Pipe_Core40",mol=mol,vmin=0,vmax=10,file_extension=file_extension)
                gauss_fitter(region="Serpens_Aquila",mol=mol,vmin=0,vmax=12,file_extension=file_extension)
                gauss_fitter(region="Serpens_MWC297",mol=mol,vmin=0,vmax=12,file_extension=file_extension)
                

def gauss_fitter(region = 'Cepheus_L1251', snr_min = 3.0, mol = 'C2S', vmin = 5.0, vmax=10.0, convolve=False, use_old_conv=False, multicore = 1, file_extension = None):

    	"""
    	Fit a Gaussian to non-NH3 emission lines from GAS.
    	It creates a cube for the best-fit Gaussian, a cube 
    	for the best-fit Gaussian with noise added back into 
    	the spectrum, and a parameter map of Tpeak, Vlsr, and FWHM
    
    	Parameters
    	----------
    	region : str
        	Name of region to reduce
    	snr_min : float
        	Lowest signal-to-noise pixels to include in the line-fitting
    	mol : str
        	name of molecule to fit
   	vmin : numpy.float
        	Minimum centroid velocity, in km/s.
    	vmax : numpy.float
        	Maximum centroid velocity, in km/s.
    	convolve : bool or float
        	If not False, specifies the beam-size to convolve the original map with
		Beam-size must be given in arcseconds
    	use_old_conv : bool
        	If True, use an already convolved map with name:
		region + '_' + mol + file_extension + '_conv.fits'
		This convolved map must be in units of km/s
    	multicore : int
		Maximum number of simultaneous processes desired
	file_extension: str
		filename extension 
    	"""
    	if file_extension:
        	root = file_extension
    	else:
        	# root = 'base{0}'.format(blorder)
        	root = 'all'

	molecules = ['C2S', 'HC7N_22_21', 'HC7N_21_20', 'HC5N']

    	MolFile = '{0}/{0}_{2}_{1}.fits'.format(region,root,mol)
	ConvFile = '{0}/{0}_{2}_{1}_conv.fits'.format(region,root,mol)
	GaussOut = '{0}/{0}_{2}_{1}_gauss_cube.fits'.format(region,root,mol)
	GaussNoiseOut = '{0}/{0}_{2}_{1}_gauss_cube_noise.fits'.format(region,root,mol)
	ParamOut = '{0}/{0}_{2}_{1}_param_cube.fits'.format(region,root,mol)

	# Load the spectral cube and convert to velocity units
	cube = SpectralCube.read(MolFile)
	cube_km = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio')

	# If desired, convolve map with larger beam 
	# or load previously created convolved cube
	if convolve:
		cube = SpectralCube.read(MolFile)
		cube_km_1 = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio')
		beam = radio_beam.Beam(major=convolve*u.arcsec, minor=convolve*u.arcsec, pa=0*u.deg)
		cube_km = cube_km_1.convolve_to(beam)
		cube_km.write(ConvFile, format='fits', overwrite=True)
	if use_old_conv:
		cube_km = SpectralCube.read(ConvFile)
	
	# Define the spectral axis in km/s
	spectra_x_axis_kms = np.array(cube_km.spectral_axis) 

	# Find the channel range corresponding to vmin and vmax
	# -- This is a hold-over from when I originally set up the code to 
	#    use a channel range rather than velocity range.
	#    Can change later, but this should work for now. 
	low_channel = np.where(spectra_x_axis_kms<=vmax)[0][0]+1 # Add ones to change index to channel
	high_channel = np.where(spectra_x_axis_kms>=vmin)[0][-1]+1 # Again, hold-over from older setup
	peak_channels = [low_channel, high_channel]

	# Create cubes for storing the fitted Gaussian profiles
	# and the Gaussians with noise added back into the spectrum
	header = cube_km.header
	cube_gauss = np.array(cube_km.unmasked_data[:,:,:])
	cube_gauss_noise = np.array(cube_km.unmasked_data[:,:,:])
	shape = np.shape(cube_gauss)

	# Set up a cube for storing fitted parameters
	param_cube = np.zeros((6, shape[1], shape[2]))
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
	# with SNR above snr_min
	x = []
	y = []
	pixels = 0
	for (i,j), value in np.ndenumerate(cube_gauss[0]):
     		spectra=np.array(cube_km.unmasked_data[:,i,j])
		if (False in np.isnan(spectra)):
     			rms = np.nanstd(np.append(spectra[0:(peak_channels[0]-1)], spectra[(peak_channels[1]+1):len(spectra)]))
     			if (max(spectra[peak_channels[0]:peak_channels[1]]) / rms) > snr_min:
            			pixels+=1
	    			x.append(i)
	    			y.append(j)
     		else:	
	    		cube_gauss[:,i,j]=nan_array
	    		param_cube[:,i,j]=nan_array2
			cube_gauss_noise[:,i,j]=nan_array
	print str(pixels) + ' Pixels above SNR=' + str(snr_min) 

	# Define a Gaussian fitting function for each pixel
	# i, j are the x,y coordinates of the pixel being fit
	def pix_fit(i,j):
		spectra = np.array(cube_km.unmasked_data[:,i,j])
		# Use the peak brightness Temp within specified channel 
		# range as the initial guess for Gaussian height
		max_ch = np.argmax(spectra[peak_channels[0]:peak_channels[1]])
		Tpeak = spectra[peak_channels[0]:peak_channels[1]][max_ch]
		# Use the velocity of the brightness Temp peak as 
		# initial guess for Gaussian mean
		vpeak = spectra_x_axis_kms[peak_channels[0]:peak_channels[1]][max_ch]
		rms = np.std(np.append(spectra[0:(peak_channels[0]-1)], spectra[(peak_channels[1]+1):len(spectra)]))
		err1 = np.zeros(shape[0])+rms
		# Create a noise spectrum based on rms of off-line channels
		# This will be added to best-fit Gaussian to obtain a noisy Gaussian 
		noise=np.random.normal(0.,rms,len(spectra_x_axis_kms))
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
			#if max(gauss)/rms < snr_min:
			#	noisy_gauss = nan_array
			#	gauss = nan_array
			#	params = nan_array2
		
		except RuntimeError:
			noisy_gauss = nan_array
			gauss = nan_array
			params = nan_array2
	
		return i, j, gauss, params, noisy_gauss

	# Parallel computation:
	nproc = multicore  # maximum number of simultaneous processes desired
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
		print str(counter) + ' of ' + str(pixels) + ' pixels completed \r',
		sys.stdout.flush()
	print "\n %f s for parallel computation." % (time.time() - tic)

	# Save final cubes
	# These will be in km/s units. 
	# Spectra will have larger values to the left, lower values to right
	cube_final_gauss = SpectralCube(data=cube_gauss, wcs=cube_km.wcs, header=cube_km.header)
	cube_final_gauss.write(GaussOut, format='fits', overwrite=True)
	cube_final_gauss_noise = SpectralCube(data=cube_gauss_noise, wcs=cube_km.wcs, header=cube_km.header)
	cube_final_gauss_noise.write(GaussNoiseOut, format='fits', overwrite=True)

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

	fits.writeto(ParamOut, param_cube, header=param_header, clobber=True)


def mask_gauss_fits(region='HC2',file_extension='all_rebase3',threshold=0.0125):
        # Mask Gaussian fit results based on some S/N and parameter checks
        if file_extension:
                root = file_extension
        else:
                root = 'all'
        line_list = ['C2S','HC5N','HC7N_22_21','HC7N_21_20','NH3_33']
        for line in line_list:
                gparamfits = '{0}/{0}_{2}_{1}_param_cube.fits'.format(region,root,line)
                gcubefits  = '{0}/{0}_{2}_{1}_gauss_cube.fits'.format(region,root,line)
                rmsfits    = '{0}/{0}_{2}_{1}_rms_QA.fits'.format(region,root,line)
                mom0fits   = '{0}/{0}_{2}_{1}_mom0_QA.fits'.format(region,root,line)
                m0sigfits   = '{0}/{0}_{2}_{1}_mom0_sigma_QA.fits'.format(region,root,line)
                maskedfits = '{0}/{0}_{2}_{1}_param_cube_masked.fits'.format(region,root,line)
                # Load pyspeckit cube and parameter cube if exists:
                if os.path.isfile(gparamfits):
                        gauss_cube = pyspeckit.Cube(gcubefits)
                        param_cube = SpectralCube.read(gparamfits)
                        param_pycube = pyspeckit.Cube(cube=param_cube)
                        rms = fits.getdata(rmsfits)
                        mom0 = fits.getdata(mom0fits)
                        m0sig = fits.getdata(m0sigfits)
                        # define mask 
                        mask3d = gauss_cube.cube > threshold
                        im_mask=np.sum(mask3d, axis=0)
                        # Defining mask for poor fits. Using moment map and moment map sigma here.
                        # Need to test on a number of regions. 
                        for ii in np.arange(im_mask.shape[1]):
                                for jj in np.arange( im_mask.shape[0]):
                                        # Make sure unobserved regions are masked
                                        if np.isnan(mom0[jj,ii]):
                                                param_pycube.cube[:,jj,ii] = np.nan
                                        else:
                                                if ((im_mask[jj,ii] == 0) or 
                                                    (param_pycube.cube[2,jj,ii] < 3*param_pycube.cube[5,jj,ii]) or
                                                    (param_pycube.cube[2,jj,ii] == 0) or (param_pycube.cube[4,jj,ii] > 0.3) or 
                                                    (mom0[jj,ii] < (4*m0sig[jj,ii]))):
                                                        param_pycube.cube[:,jj,ii] = 0
                        # create masked parameter cube
                        masked_param_hdu = fits.PrimaryHDU(param_pycube.cube,header=param_cube.header)
                        masked_param_hdu.writeto(maskedfits,overwrite=True)


def mask_gauss_fits_all(file_extension='all_rebase3',threshold=0.0125):
        # Get list of regions - run from images/ directory
        # Assume directories correspond to regions to be imaged
        # Update - use catalog?
        # Note: need to have already updated moment maps based on fits for this step because we're using the mom0_QA files...
        region_list = glob.glob("*/")
        for i in range(len(region_list)):
                region_list[i] = region_list[i].strip("/")
        if 'figures' in region_list: region_list.remove('figures')

        for region in region_list:
                mask_gauss_fits(region=region,file_extension=file_extension,threshold=threshold)


### Examples ###
# Fit the HC5N data in Cepheus_L1251, without convolution
#gauss_fitter(region = 'Cepheus_L1251', snr_min = 7.0, mol = 'HC5N', vmin=-6.3, vmax=-2.2, multicore=3)

# Convolve the HC5N data in Cepheus_L1251 to a spatial resolution of 64 arcseconds,
# then fit a Gaussian to all pixels above SNR=3
#gauss_fitter(region = 'Cepheus_L1251', direct = '/Users/jkeown/Desktop/GAS_dendro/', snr_min = 3.0, mol = 'HC5N', peak_channels = [402,460], convolve=64., use_old_conv=False)


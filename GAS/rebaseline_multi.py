import numpy
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares as lsq
import numpy.polynomial.legendre as legendre
import pprocess
import time
import sys
import skimage
from keystone.baseline import mad1d, legendreLoss

def get_mask(spectra, mask_percent=0.4, window_size=31):
      	"""  
 Returns a mask of channels to be used for baseline fitting.
 Function calculates the standard deviation of a 31 pixel window
 centred on each pixel. A percentage of the pixels with the lowest 
 standard deviation for their window are then chosen for baseline fitting.    
   
 spectra = input spectra as numpy array
 mask_percent = percentage of pixels to select for baseline fitting
 window_size = width of window, in units of channels, centred on each channel    
       	""" 
	spec_len = len(spectra)
	sample = int(spec_len*mask_percent)
	left = np.zeros(int(window_size/2))+np.nanstd(spectra[0:window_size]) # For the first few channels without a full window, use first window
	right = np.zeros(int(window_size/2))+np.nanstd(spectra[-window_size:]) # For the last few channels without a full window, use last window
	middle = np.nanstd(skimage.util.view_as_windows(spectra, window_size, 1),axis=1)
	stds = np.concatenate((left, middle, right))
	mask = np.arange(spec_len)[np.argsort(stds)[:sample]]
	#median_std = np.median(stds)*3
	#mask = np.where(np.array(stds)<median_std)[0]
	#mask = np.concatenate((mask, np.arange(-50, 50)))
	#if 1>0: #len(mask)<30
		#print 'yes'
		#mask = np.arange(len(spectra))[np.argsort(stds)[:500]]
		#mask=np.arange(len(spectra))[0::5]

	#plt.plot(range(len(spectra)), spectra)
	#plt.plot(np.arange(len(spectra))[mask], spectra[mask])
	#plt.show()
	return mask

def redchisqg(ydata,ymod,deg=2,sd=None):     
      """Returns the reduced chi-square error statistic for an arbitrary model,   
 chisq/nu, where nu is the number of degrees of freedom. If individual   
 standard deviations (array sd) are supplied, then the chi-square error   
 statistic is computed as the sum of squared errors divided by the standard   
 deviations. See http://en.wikipedia.org/wiki/Goodness_of_fit for reference.  
   
 ydata,ymod,sd assumed to be Numpy arrays. deg integer.  
   
 Usage:  
 >>> chisq=redchisqg(ydata,ymod,n,sd)  
 where  
  ydata : data  
  ymod : model evaluated at the same x points as ydata  
  n : number of free parameters in the model  
  sd : uncertainties in ydata  
   
 Rodrigo Nemmen  
 http://goo.gl/8S1Oo"""  
      # Chi-square statistic  
      if sd==None:  
           chisq=np.sum((ydata-ymod)**2)  
      else:  
           chisq=np.sum( ((ydata-ymod)/sd)**2 )  
             
      # Number of degrees of freedom assuming 2 free parameters  
      nu=ydata.size-1-deg  
        
      return chisq/nu

def get_chi(blorder_max, ydata, xdata, blindex, noise):
    """
 Returns the best-fit Legendre polynomial to an input spectra, 
 based on a comparison of the reduced chi-squared values for each order 
 of the polynomial.     
   
 blorder_max = maximum order of polynomial to fit
 ydata = spectrum y values
 xdata = spectrum x values
 blindex = the indices of the spectra to include in baseline fitting
 noise = rms noise of spectrum"""
    opts = lsq(legendreLoss, np.zeros(blorder_max + 1), args=(ydata[blindex],
                                                          xdata[blindex],
                                                        noise), loss='arctan')
    chis = [] # first entry is order=1, second is order=2, ..., up to order=blorder
    ymods = []
    # Loop through each order of polynomial and create model baseline
    # and calculate that model's chi-squared values
    for i in np.arange(blorder_max)+1:
	# i is the polynomial order
    	ymod = legendre.legval(xdata, np.array(opts.x)[0:i+1]) #+1 to get last term
	ymods.append(ymod)
    	chi = redchisqg(ydata,ymod,deg=i+1)
	chis.append(chi)
    # Find the model with the lowest chi-squared value
    find_low = np.where(np.array(chis)==min(chis))
    low_model = np.array(ymods)[find_low]
    return ymod

def robustBaseline_chi(y, baselineIndex, blorder_max=3, noiserms=None):
    """  
 Returns a baseline subtracted spectrum, based on the best-fitting polynomial 
 for a range of polynomial orders.    
   
 y = input spectra
 baselineIndex = indices of spectra to include in baseline fitting
 blorder_max = largest order polynomial to fit (fit from blorder_max down to order of 1) 
 noiserms = rms noise of spectrum"""
    x = np.linspace(-1, 1, len(y))
    if noiserms is None:
        noiserms = mad1d((y - np.roll(y, -2))[baselineIndex]) * 2**(-0.5)
    low_model = get_chi(blorder_max=blorder_max, ydata=y, xdata=x, blindex=baselineIndex, noise=noiserms)
    #plt.plot(range(len(y)), y)
    #plt.plot(np.arange(len(y))[baselineIndex], y[baselineIndex])  
    #plt.plot(range(len(x)), low_model, color='red')
    #plt.show()
    return y - low_model

def rebase(ch, i, data, mask_percent=0.4, blorder_max=3, window_size=31):
	"""  
 Parallelizable function to feed into pprocess. Returns a baseline-subtracted
 spectrum and its indices on the image plane.    
   
 i,j = x and y indices for spectrum location on image plane
 spectra = spectrum on which a baseline will be subtracted
 mask_percent = percentage of pixels to select for baseline fitting
 blorder_max = largest order polynomial to fit (fit from blorder_max down to order of 1) 
	"""
	for ii in i:
		for j in range(data.shape[2]):
			spectra=np.array(data[:,ii,j])
			if (False in np.isnan(spectra)): #and (m/std > 10.):
				mask = get_mask(spectra, mask_percent=mask_percent, window_size=window_size)
				spectra = robustBaseline_chi(spectra, mask, blorder_max=blorder_max, noiserms=None)
			ch.send((ii, j, spectra))

def rebase_multi(filename, nproc=8, mask_percent=0.4, blorder_max=3, window_size=31):
	"""  
 Returns a baseline-subtracted cube. Can be run with parallel processes.    
   
 filename = name of datacube to process (including its path)
 nproc = number of parallel processes desired
 mask_percent = percentage of pixels to select for baseline fitting
 blorder_max = largest order polynomial to fit (fit from blorder_max down to order of 1) 
	"""
	cube = SpectralCube.read(filename)

	queue = pprocess.Queue(limit=nproc, continuous=1)
	calc = queue.manage(rebase)
	tic = time.time()

	# create cube to store rebaselined data
	cube_out = np.zeros(cube.shape) * np.nan
	pixels = cube.shape[1] * cube.shape[2]

	counter = 0
	for i in np.array_split(range(cube.shape[1]), nproc):
		calc(i, data=cube, mask_percent=mask_percent, blorder_max=blorder_max, window_size=window_size)

	for i, j, ss in queue:
		cube_out[:,i,j]=ss
		counter+=1
		print str(counter) + ' of ' + str(pixels) + ' pixels completed \r',
		sys.stdout.flush()
	print "\n %f s for parallel computation." % (time.time() - tic)
	
	cube_final = SpectralCube(data=cube_out, wcs=cube.wcs, header=cube.header)
	cube_final.write(filename[0:-5] + '_rebase_multi.fits', format='fits', overwrite=True)

# Usage examples:
#rebase_multi('/lustre/pipeline/scratch/KEYSTONE/images/W48_NH3_11_sess31-sess31.fits', nproc=4)
#rebase_multi('/lustre/pipeline/scratch/KEYSTONE/images/NGC2264_NH3_11_all.fits', mask_percent=0.2)
#rebase_multi('/lustre/pipeline/scratch/KEYSTONE/images/Rosette_NH3_11_all.fits', mask_percent=0.4)
#rebase_multi('/lustre/pipeline/scratch/KEYSTONE/images/NGC2264_NH3_22_all.fits', mask_percent=0.2, blorder_max=1)

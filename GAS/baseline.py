import numpy as np
import numpy.polynomial.legendre as legendre
import warnings
from .utils import VlsrByCoord
import catalogs
from scipy.optimize import least_squares as lsq
from spectral_cube import SpectralCube
import astropy.units as u
import itertools
import astropy.utils.console as console
import pyspeckit.spectrum.models.ammonia_constants as acons
#from pyspeckit.parallel_map import parallel_map

def ammoniaWindow(spectrum, spaxis, freqthrow = 4.11*u.MHz, 
                  window = 2, v0 = 8.5, line='oneone',outerwindow=None):
    mask = np.zeros_like(spectrum,dtype = np.bool)
    voffs = np.array([dv for dv in acons.voff_lines_dict[line]])
    for hyperfine in voffs:
        mask[(spaxis>(v0+voff-window))*(spaxis<(v0+voff+window))] = True
    deltachan = freqthrow/((spaxis[1]-spaxis[0])/299792.458*acons.freq_dict[line]*u.Hz)
    deltachan = deltachan.to(u.dimensionless_unscaled).value
    deltachan = (np.floor(np.abs(
                deltachan.to(u.dimensionless_unscaled).value))).astype(np.int)
    mask = np.logical_or(mask,np.r_[mask[deltachan:-1],
                                    np.zeros(deltachan+1,dtype=np.bool)])
    mask = np.logical_or(mask,np.r_[np.zeros(deltachan+1,dtype=np.bool),
                                    mask[0:(-deltachan-1)]])
    # TODO get in frequency switch throw
    if outerwindow is not None:
        mask[(spaxis>(v0+outerwindow+voffs.max()))] = True
        mask[(spaxis<(v0-outerwindow-voffs.min()))] = True

    return ~mask

def tightWindow(spectrum, spaxis,
                window = 5,
                outerwindow = None,
                v0 = 8.5,freqthrow = 4.11*u.MHz):
    mask = np.zeros_like(spectrum,dtype = np.bool)
    mask[(spaxis>(v0-window))*(spaxis<(v0+window))]=True
    deltachan = freqthrow/((spaxis[1]-spaxis[0])/299792.458*23.5*u.GHz)
    deltachan = (np.floor(np.abs(
                deltachan.to(u.dimensionless_unscaled).value))).astype(np.int)
    mask = np.logical_or(mask,np.r_[mask[deltachan:-1],\
                                        np.zeros(deltachan+1,dtype=np.bool)])
    mask = np.logical_or(mask,np.r_[np.zeros(deltachan+1,dtype=np.bool),
                                    mask[0:(-deltachan-1)]])
    if outerwindow is not None:
        mask[(spaxis>(v0+outerwindow))] = True
        mask[(spaxis<(v0-outerwindow))] = True
    return(~mask)

def mad1d(x):
    med0 = np.median(x)
    return np.median(np.abs(x-med0))*1.4826

def legendreLoss(coeffs,y,x,noise):
    return (y-legendre.legval(x,coeffs))/noise

def robustBaseline(y,baselineIndex,blorder=1,noiserms=None):
    x = np.linspace(-1,1,len(y))
    if noiserms is None:
        noiserms = mad1d((y-np.roll(y,-2))[baselineIndex])
    opts = lsq(legendreLoss,np.zeros(blorder+1),args=(y[baselineIndex],
                                                      x[baselineIndex],
                                                      noiserms),loss='arctan')
#    import pdb; pdb.set_trace()

    return y-legendre.legval(x,opts.x)

def rebaseline(filename, blorder = 1, 
               baselineRegion = [slice(0,262,1),slice(-512,0,1)], 
               windowFunction = None, **kwargs):
    """
    Rebaseline a data cube using robust regression of Legendre polynomials.

    Parameters
    ----------
    filename : string
        FITS filename of the data cube
    blorder : int
        Order of the polynomial to fit to the data
    baselineRegion : list
        List of slices defining the default region of the spectrum, in 
        channels, to be used for the baseline fitting. 
    windowFuntion : function
        Name of function to be used that will accept spectrum data, and velocity
        axis and will return a binary mask of the channels to be used in the 
        baseline fitting.  Extra **kwargs are passed to windowFunction to do with
        as it must.

    Returns
    -------
    Nothing.  A new FITS file is written out with the suffix 'rebaseN' where N 
    is the baseline order

    """
    cube = SpectralCube.read(filename)
    originalUnit = cube.spectral_axis.unit
    cube = cube.with_spectral_unit(u.km/u.s,velocity_convention='radio')
    
    goodposition = np.isfinite(cube.apply_numpy_function(np.max,axis=0))
    y,x = np.where(goodposition)
    outcube = np.zeros(cube.shape)*np.nan
    flatshape = (cube.shape[0],cube.shape[1]*cube.shape[2])

    RegionName = (filename.split('_'))[0]

    if hasattr(windowFunction,'__call__'):
        catalog = catalogs.GenerateRegions()

    nuindex = np.arange(cube.shape[0])
    runmin = nuindex[-1]
    runmax = nuindex[0]
    for thisy,thisx in console.ProgressBar(zip(y,x)):
        spectrum = cube[:,thisy,thisx].value

        if hasattr(windowFunction,'__call__'):
            _, Dec, RA = cube.world[0,thisy,thisx]
            # This determines a v0 appropriate for the region
            v0 = VlsrByCoord(RA.value,Dec.value,RegionName,regionCatalog = catalog)
            baselineIndex = windowFunction(spectrum,
                                           cube.spectral_axis.to(u.km/u.s).value,
                                           v0 = v0, **kwargs)
        else:
            baselineIndex = np.concatenate([nuindex[ss] for ss in baselineRegion])
        runmin = np.min([nuindex[baselineIndex].min(),runmin])
        runmax = np.max([nuindex[baselineIndex].max(),runmax])
        # Use channel-to-channel difference as the noise value.
        noise = mad1d((spectrum-np.roll(spectrum,-2))[baselineIndex])*2**(-0.5)

        outcube[:,thisy,thisx] = robustBaseline(spectrum,baselineIndex,
                                                blorder=blorder,noiserms=noise)
    outsc = SpectralCube(outcube,cube.wcs,header=cube.header)
    outsc = outsc[runmin:runmax,:,:] # cut beyond baseline edges
    # Return to original spectral unit
    outsc = outsc.with_spectral_unit(originalUnit)
    outsc.write(filename.replace('.fits','_rebase{0}.fits'.format(blorder)),
                overwrite=True)
    

import numpy as np
import numpy.polynomial.legendre as legendre
import warnings
from .intervals import VelocitySet,VelocityInterval
from .utils import VlsrByCoord
import catalogs
from scipy.optimize import least_squares as lsq
from spectral_cube import SpectralCube
import astropy.units as u
import itertools
import astropy.utils.console as console
import pyspeckit.spectrum.models.ammonia_constants as acons
#from pyspeckit.parallel_map import parallel_map

oneone = [VelocityInterval(voff-0.001,voff+0.001) for 
          voff in acons.voff_lines_dict['oneone']]
twotwo = [VelocityInterval(voff-0.001,voff+0.001) for 
          voff in acons.voff_lines_dict['twotwo']]
threethree = [VelocityInterval(voff-0.001,voff+0.001) for 
              voff in acons.voff_lines_dict['threethree']]

def ammoniaWindow(spectrum, spaxis, freqthrow = 4.11*u.MHz, 
                  window = 2, v0 = 8.5, line='oneone'):
    mask = np.zeros_like(spectrum,dtype = np.bool)
    for hyperfine in acons.voff_lines_dict[line]:
        mask[(spaxis>(v0+voff-window))*(spaxis<(v0+voff+window))] = True
    deltachan = freqthrow/((spaxis[1]-spaxis[0])/299792.458*acons.freq_dict[line])
    # TODO get in frequency switch throw
    return ~mask

def tightWindow(spectrum, spaxis,
                window = 5,
                v0 = 8.5):
    mask = np.zeros_like(spectrum,dtype = np.bool)
    mask[(spaxis>(v0-window))*(spaxis<(v0+window))]=True
    return(~mask)

def mad1d(x):
    med0 = np.median(x)
    return np.median(np.abs(x-med0))*1.4826

def legendreLoss(coeffs,y,x,noise):
    return np.sum(np.abs((y-legendre.legval(x,coeffs))/noise))

def robustBaseline(y,baselineIndex,blorder=1):
    x = np.linspace(-1,1,len(y))
    noise = mad1d((y-np.roll(y,-2))[baselineIndex])*2**(-0.5)
    opts = lsq(legendreLoss,np.zeros(blorder+1),args=(y[baselineIndex],
                                                      x[baselineIndex],
                                                      noise),loss='arctan')
    return y-legendre.legval(x,opts.x)

def rebaseline(filename, blorder = 1, 
               baselineRegion = [slice(0,262,1),slice(-512,0,1)], 
               windowFunction = None, **kwargs):

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
            v0 = VlsrByCoord(RA.value,Dec.value,RegionName,regionCatalog = catalog)
            baselineIndex = windowFunction(spectrum,
                                           cube.spectral_axis.to(u.km/u.s).value,
                                           v0 = v0, **kwargs)
        else:
            baselineIndex = np.concatenate([nuindex[ss] for ss in baselineRegion])
        runmin = np.min([nuindex[baselineIndex].min(),runmin])
        runmax = np.max([nuindex[baselineIndex].max(),runmax])

        outcube[:,thisy,thisx] = robustBaseline(spectrum,baselineIndex,blorder=blorder)
#        import pdb; pdb.set_trace()
    outsc = SpectralCube(outcube,cube.wcs,header=cube.header)
#    outsc = outsc[runmin:runmax,:,:] # cut baseline edges
    outsc.with_spectral_unit(originalUnit)
    outsc.write(filename.replace('.fits','_rebase{0}.fits'.format(blorder)),overwrite=True)
    

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

def tightWindow(spectrum, spaxis,
                vwindow = 5,
                v0 = 8.5):

    EmissionWindow = VelocityInterval(-vwindow,vwindow)
    EmissionWindow.applyshift(v0)
    BaselineWindow = VelocitySet([VelocityInterval(spaxis.min(),spaxis.max())])-\
        EmissionWindow
    slices = BaselineWindow.toslice(spaxis)
    return(slices)

def mad1d(x):
    med0 = np.median(x)
    return np.median(np.abs(x-med0))*1.4826

def legendreLoss(coeffs,y,x,noise):
    return np.sum(np.abs((y-legendre.legval(x,coeffs))/noise))

def robustBaseline(y,baselineIndex,blorder=1):
    x = np.linspace(-1,1,len(y))
    noise = mad1d((y-np.roll(y,-2))[baselineIndex])*2**(-0.5)
    opts = lsq(legendreLoss,np.zeros(blorder+1),args=(y,x,noise),loss='arctan')
    return y-legendre.legval(x,opts.x)

def rebaseline(filename, blorder = 1, 
               baselineRegion = [slice(0,262,1),slice(-512,0,1)], 
               baselineFunction = None, **kwargs):

    cube = SpectralCube.read(filename)
    cube = cube.with_spectral_unit(u.km/u.s,velocity_convention='radio')

    goodposition = np.isfinite(cube.apply_numpy_function(np.max,axis=0))
    y,x = np.where(goodposition)
    outcube = np.zeros(cube.shape)*np.nan

    RegionName = (filename.split('_'))[0]

    if hasattr(baselineFunction,'__call__'):
        catalog = catalogs.GenerateRegions()

    for thisy,thisx in itertools.izip(y,x):
        spectrum = cube[:,thisy,thisx].value
        nuindex = np.arange(len(spectrum))
        if hasattr(baselineFunction,'__call__'):
            _, Dec, RA = cube.world[0,thisy,thisx]
            v0 = VlsrByCoord(RA.value, Dec.value,RegionName,regionCatalog = catalog)
            baselineRegion = baselineFunction(spectrum,
                                              cube.spectral_axis.to(u.km/u.s).value,
                                              v0 = v0, **kwargs)
        baselineIndex = np.concatenate([nuindex[ss] for ss in baselineRegion])
        outcube[:,thisy,thisx] = robustBaseline(spectrum,baselineIndex,blorder=blorder)
        
    outsc = SpectralCube(outcube,cube.wcs,header=cube.header)
    outsc = outsc[262:-512,:,:] # cut baseline edges
    outsc.write(filename.replace('.fits','_rebase{0}.fits'.format(blorder)),overwrite=True)
    
    
    # EW = VelocityInterval(-5e3,5e3)
    # EW.applyshift(8.500e3)
    # if hasattr(baselineGenerator,'__call__'):
    #     baselineRegion = baselineGenerator(spectrum,EmissionWindow = EW)
    #     baselineIndex = np.concatenate([nuindex[ss] for ss in baselineRegion])
    # else:

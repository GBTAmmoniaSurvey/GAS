import numpy as np
import numpy.polynomial.legendre as legendre
import warnings
from .intervals import VelocitySet,VelocityInterval
from scipy.optimize import least_squares as lsq
from spectral_cube import SpectralCube
import itertools

def testBLgen(spectrum,
              EmissionWindow = VelocityInterval(-5e3,5e3),
              v0 = 8.5e3):
    """
    Returns good slices for baseline fitting given an input GBTIDL FITS structure
    """
    
    cms = 299792458.
    # First calculate full velocity range of spectra in LSRK frame
    nChannels = len(spectrum['DATA'])
    lowedge = (1-(spectrum['CRVAL1']+spectrum['CDELT1']*(1-spectrum['CRPIX1']))/
               spectrum['RESTFREQ'])*cms-spectrum['VFRAME']
    hiedge = (1-(spectrum['CRVAL1']+spectrum['CDELT1']*(nChannels-spectrum['CRPIX1']))/
              spectrum['RESTFREQ'])*cms-spectrum['VFRAME']
    lowedge +=20000
    hiedge -= 20000
    
    EmissionWindow = VelocitySet([VelocityInterval(-70e3,70e3)])
#    v0 = 8500
#    EmissionWindow.applyshift(v0)
    BaselineWindow = VelocitySet([VelocityInterval(lowedge,hiedge)])-\
        EmissionWindow
    slices = BaselineWindow.toslice(cdelt = spectrum['CDELT1'], 
                                    crval = spectrum['CRVAL1'],
                                    vframe = spectrum['VFRAME'], 
                                    crpix = spectrum['CRPIX1'],
                                    restfreq = spectrum['RESTFREQ'])
    return(slices)

def legendreLoss(y,x,noise,coeffs):
    return np.abs((y-legendre.legval(x,coeffs))/noise)

def robustBaseline(y,baselineIndex):
    x = np.linspace(-1,1,len(y))
    



def rebaseline(filename, blorder = 1, baselineGenerator = None,        
               baselineRegion = [slice(0,262,1),slice(-512,0,1)], baselineFunction = ):


    cube = SpectralCube.read(filename)
    goodposition = cube.apply_numpy_function(np.any,axis=0)
    y,x = np.where(goodposition)
    for thisy,thisx in itertools.izip(y,x):
        spectrum = cube[:,thisy,thisx]
        


    
    
    # EW = VelocityInterval(-5e3,5e3)
    # EW.applyshift(8.500e3)
    # if hasattr(baselineGenerator,'__call__'):
    #     baselineRegion = baselineGenerator(spectrum,EmissionWindow = EW)
    #     baselineIndex = np.concatenate([nuindex[ss] for ss in baselineRegion])
    # else:

#from sdpy import makecube
import numpy as np
import glob
from astropy.io import fits
import astropy.wcs as wcs
import itertools
from scipy.special import j1
import pdb
import numpy.fft as fft
import astropy.utils.console as console
import astropy.units as u
import astropy.constants as con
import numpy.polynomial.legendre as legendre
import warnings
# import baseline
import gbtpipe
from . import __version__
from .utils import VlsrByCoord

def baselineSpectrum(spectrum, order=1, baselineIndex=()):
    x = np.linspace(-1, 1, len(spectrum))
    coeffs = legendre.legfit(x[baselineIndex], spectrum[baselineIndex], order)
    spectrum -= legendre.legval(x, coeffs)
    return(spectrum)


def freqShiftValue(freqIn, vshift, convention='RADIO'):
    cms = 299792458.
    if convention.upper() in 'OPTICAL':
        return freqIn / (1.0 + vshift / cms)
    if convention.upper() in 'TRUE':
        return freqIn * ((cms + vshift) / (cms - vshift))**0.5
    if convention.upper() in 'RADIO':
        return freqIn * (1.0 - vshift / cms)


def channelShift(x, ChanShift):
    # Shift a spectrum by a set number of channels.
    ftx = np.fft.fft(x)
    m = np.fft.fftfreq(len(x))
    phase = np.exp(2 * np.pi * m * 1j * ChanShift)
    x2 = np.real(np.fft.ifft(ftx * phase))
    return(x2)


def jincGrid(xpix, ypix, xdata, ydata, pixPerBeam=None):
    a = 1.55 / (3.0 / pixPerBeam)
    b = 2.52 / (3.0 / pixPerBeam)

    Rsup = 1.0 * pixPerBeam  # Support radius is 1 FWHM
    dmin = 1e-4
    dx = (xdata - xpix)
    dy = (ydata - ypix)

    pia = np.pi / a
    b2 = 1. / (b**2)
    distance = np.sqrt(dx**2 + dy**2)

    ind  = (np.where(distance <= Rsup))
    d = distance[ind]
    wt = j1(d * pia) / \
        (d * pia) * \
        np.exp(-d**2 * b2)
#    wt[ind] = np.exp(-distance[ind]**2*b2)*\
#              np.sin(pia*distance[ind])/\
#              (pia*distance[ind])
    wt[(d < dmin)] = 0.5  # Peak of the jinc function is 0.5 not 1.0
    return(wt, ind)


def autoHeader(filelist, beamSize=0.0087, pixPerBeam=3.0):
    RAlist = []
    DEClist = []
    for thisfile in filelist:
        s = fits.getdata(thisfile)
        try:
            RAlist = RAlist + [s['CRVAL2']]
            DEClist = DEClist + [s['CRVAL3']]
        except:
            pdb.set_trace()

    longitude = np.array(list(itertools.chain(*RAlist)))
    latitude = np.array(list(itertools.chain(*DEClist)))
    longitude = longitude[longitude != 0]
    latitude = latitude[latitude != 0]
    minLon = longitude.min()
    maxLon = longitude.max()
    minLat = latitude.min()
    maxLat = latitude.max()

    naxis2 = np.ceil((maxLat - minLat) /
                     (beamSize / pixPerBeam) + 2 * pixPerBeam)
    crpix2 = naxis2 / 2
    cdelt2 = beamSize / pixPerBeam
    crval2 = (maxLat + minLat) / 2
    ctype2 = s[0]['CTYPE3'] + '--TAN'
    # Negative to go in the usual direction on sky:
    cdelt1 = -beamSize / pixPerBeam

    naxis1 = np.ceil((maxLon - minLon) /
                     (beamSize / pixPerBeam) *
                     np.cos(crval2 / 180 * np.pi) + 2 * pixPerBeam)
    crpix1 = naxis1 / 2
    crval1 = (minLon + maxLon) / 2
    ctype1 = s[0]['CTYPE2'] + '---TAN'
    outdict = {'CRVAL1': crval1, 'CRPIX1': crpix1,
               'CDELT1': cdelt1, 'NAXIS1': naxis1,
               'CTYPE1': ctype1, 'CRVAL2': crval2,
               'CRPIX2': crpix2, 'CDELT2': cdelt2,
               'NAXIS2': naxis2, 'CTYPE2': ctype2}

    return(outdict)


def addHeader_nonStd(hdr, beamSize, Data_Unit):
    if Data_Unit == 'Tmb':
        hdr['BUNIT'] = 'K'
    hdr['BMAJ'] = beamSize
    hdr['BMIN'] = beamSize
    hdr['BPA'] = 0.0
    hdr['TELESCOP'] = 'GBT'
    hdr['INSTRUME'] = 'KFPA'
    return(hdr)


def griddata(templateHeader=None,
             gridFunction=jincGrid,
             rootdir='/lustre/pipeline/scratch/GAS/',
             region='NGC1333',
             dirname='NGC1333_NH3_11',
             startChannel=762, endChannel=3584,
             doBaseline=True,
             baselineRegion=[slice(762, 1024, 1), slice(3072, 3584, 1)],
             blorder=1,
             Sessions=None,
             file_extension=None,
             flagRMS=True,
             flagSpike=True,
             **kwargs):

    if not Sessions:
        filelist = glob.glob(rootdir + '/' + region +
                             '/' + dirname + '/*fits')
        if not file_extension:
            file_extension = '_all'
        history_message = 'Gridding of data using all sessions'
    else:
        filelist = []
        for scan_i in Sessions:
                filelist.extend(glob.glob(rootdir + '/' + region +
                                          '/' + dirname + '/*_sess' +
                                          str(scan_i) + '.fits'))
        if isinstance(Sessions, list):
            if not file_extension:
                file_extension = '_sess{0}-sess{1}'.format(Sessions[0],
                                                           Sessions[-1])
            if (Sessions[-1] + 1. - Sessions[0]) / len(Sessions) == 1.0:
                history_message = 'Gridding of data using sessions' \
                    'between {0} and {1}'.format(Sessions[0], Sessions[-1])
            else:
                history_message = 'Gridding of data using sessions: '
                for scan_i in Sessions:
                    history_message += '{0}, '.format(scan_i)
        else:
            if not file_extension:
                file_extension = '_sess{0}'.format(Sessions)
            history_message = 'Gridding of data using session'\
                '{0}'.format(Sessions)

    if len(filelist) == 0:
        warnings.warn('There are no FITS files to process '
                      'in ' + rootdir + '/' + region + '/' + dirname)
        return
    # check that every file in the filelist is valid
    # If not then remove it and send warning message
    for file_i in filelist:
        try:
            fits.open(file_i)
        except:
            warnings.warn('file {0} is corrupted'.format(file_i))
            filelist.remove(file_i)
    outdir= rootdir + '/images/' + region + '/' 
    outname = dirname + file_extension
    gbtpipe.Gridding.griddata(filelist,
                              startChannel=startChannel,
                              endChannel=endChannel,
                              doBaseline=doBaseline,
                              baselineRegion=baselineRegion,
                              blorder=blorder, rebaseorder=3,
                              flagRMS=flagRMS,
                              outdir=outdir, outname=outname,
                              VlsrByCoord=VlsrByCoord,
                              flagSpike=flagSpike, **kwargs)

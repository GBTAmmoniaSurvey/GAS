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
from .intervals import VelocitySet,VelocityInterval
from . import __version__

def baselineSpectrum(spectrum,order=1,baselineIndex=()):
    x=np.arange(len(spectrum))-len(spectrum)/2
    coeffs = legendre.legfit(x[baselineIndex],spectrum[baselineIndex],order)
    spectrum -= legendre.legval(x,coeffs)
    return(spectrum)

def freqShiftValue(freqIn,vshift,convention='RADIO'):
    cms = 299792458.
    if convention.upper() in 'OPTICAL':
        return freqIn/(1.0+vshift/cms)
    if convention.upper() in 'TRUE':
        return freqIn*((cms+vshift)/(cms-vshift))**0.5
    if convention.upper() in 'RADIO':
        return freqIn*(1.0-vshift/cms)
    
def channelShift(x,ChanShift):
# Shift a spectrum by a set number of channels.  
    ftx = np.fft.fft(x)
    m = np.fft.fftfreq(len(x))
    phase = np.exp(2*np.pi*m*1j*ChanShift)
    x2 = np.real(np.fft.ifft(ftx*phase))
    return(x2)

def jincGrid(xpix,ypix,xdata,ydata, pixPerBeam = None):
    a = 1.55/(3.0/pixPerBeam)
    b = 2.52/(3.0/pixPerBeam)
    
    Rsup = 1.0*pixPerBeam # Support radius is 1 FWHM
    dmin = 1e-4
    dx = (xdata - xpix)#/pixPerBeam
    dy = (ydata - ypix)#/pixPerBeam

    pia = np.pi/a
    b2 = 1./(b**2)
    distance = np.sqrt(dx**2+dy**2)
    
    ind  = (np.where(distance<=Rsup))
    d = distance[ind]
    wt = j1(d*pia)/\
        (d*pia)*\
        np.exp(-d**2*b2)
#    wt[ind] = np.exp(-distance[ind]**2*b2)*\
#              np.sin(pia*distance[ind])/\
#              (pia*distance[ind])
    wt[(d<dmin)]=0.5  #Peak of the jinc function is 0.5 not 1.0
    return(wt,ind)

def testBLgen(spectrum,
              EmissionWindow = VelocityInterval(-5e3,5e3),
              v0 = 8.5e3):
    """
    Returns good slices for baseline fitting given an input GBTIDL FITS structure
    """
    
    cms = 299792458.
    # First calculate full velocity range of spectra in LSRK frame
    nChannels = len(spectrum['DATA'])
#    lowedge = (1-(spectrum['CRVAL1']+spectrum['CDELT1']*(1-spectrum['CRPIX1']))/
#               spectrum['RESTFREQ'])*cms-spectrum['VFRAME']
#    hiedge = (1-(spectrum['CRVAL1']+spectrum['CDELT1']*(nChannels-spectrum['CRPIX1']))/
#              spectrum['RESTFREQ'])*cms-spectrum['VFRAME']
    lowedge = -100000
    hiedge = 100000
    EmissionWindow = VelocityInterval(-25e3,25e3)-VelocityInterval(-14e3,-12e3)
    EmissionWindow = EmissionWindow-VelocityInterval(12e3,14e3)
    coeffs = [-2.8256074 , -4.65791997,  9.14502305]
    v0 = coeffs[2]+\
        coeffs[0]*(spectrum['CRVAL2']-83.446122802665869)+\
        coeffs[1]*(spectrum['CRVAL3']+6.0050560818354661)
    EmissionWindow.applyshift(v0)
    BaselineWindow = VelocitySet([VelocityInterval(lowedge,hiedge)])-\
        EmissionWindow
    slices = BaselineWindow.toslice(cdelt = spectrum['CDELT1'], 
                                    crval = spectrum['CRVAL1'],
                                    vframe = spectrum['VFRAME'], 
                                    crpix = spectrum['CRPIX1'],
                                    restfreq = spectrum['RESTFREQ'])
    return(slices)


def autoHeader(filelist, beamSize = 0.0087, pixPerBeam = 3.0):
    RAlist = []
    DEClist = []
    for thisfile in filelist:
        s = fits.getdata(thisfile)
        try:
            RAlist = RAlist + [s['CRVAL2']]
            DEClist = DEClist +[s['CRVAL3']]
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

    naxis2 = np.ceil((maxLat-minLat)/(beamSize/pixPerBeam)+2*pixPerBeam)
    crpix2 = naxis2/2
    cdelt2 = beamSize/pixPerBeam
    crval2 = (maxLat+minLat)/2
    ctype2 = s[0]['CTYPE3']+'--TAN'
    cdelt1 = -beamSize/pixPerBeam # Negative to go in the usual direction on sky.
    naxis1 = np.ceil((maxLon-minLon)/(beamSize/pixPerBeam)*\
                         np.cos(crval2/180*np.pi)+2*pixPerBeam)
    crpix1 = naxis1/2
    crval1 = (minLon+maxLon)/2
    ctype1 = s[0]['CTYPE2']+'---TAN'
    outdict = {'CRVAL1':crval1,'CRPIX1':crpix1,'CDELT1':cdelt1,'NAXIS1':naxis1,'CTYPE1':ctype1,
               'CRVAL2':crval2,'CRPIX2':crpix2,'CDELT2':cdelt2,'NAXIS2':naxis2,'CTYPE2':ctype2}

    return(outdict)

def addHeader_nonStd( hdr, beamSize, Data_Unit):
    if Data_Unit == 'Tmb':
        hdr['BUNIT'] = 'K'
    hdr['BMAJ'] = beamSize
    hdr['BMIN'] = beamSize
    hdr['BPA'] = 0.0
    hdr['TELESCOP']='GBT'
    hdr['INSTRUME']='KFPA'
    return(hdr)

def griddata(pixPerBeam = 3.0,
             templateHeader = None,
             gridFunction = jincGrid, 
             rootdir = '/lustre/pipeline/scratch/GAS/',
             region = 'NGC1333',
             dirname = 'NGC1333_NH3_11',
             startChannel = 1024, endChannel = 3072,
             doBaseline = True,
             baselineRegion = [slice(512,1024,1),slice(3072,3584,1)],
             blorder = 1,
             baselineGenerator = None,
             Sessions = None,
             file_extension = None):

    if not Sessions:
        filelist = glob.glob(rootdir+'/'+region+'/'+dirname+'/*fits')
        if not file_extension:
            file_extension='_all'
        history_message='Gridding of data using all sessions'
    else:
        filelist = []
        for scan_i in Sessions:
                filelist.extend(glob.glob(rootdir+'/'+region+'/'+dirname+'/*_sess'+str(scan_i)+'.fits'))
        if isinstance(Sessions, list):
            if not file_extension:
                file_extension='_sess{0}-sess{1}'.format(Sessions[0],Sessions[-1])
            if (Sessions[-1]+1.-Sessions[0])/len(Sessions) == 1.0:
                history_message='Gridding of data using sessions between {0} and {1}'.format(Sessions[0],Sessions[-1])
            else:
                history_message='Gridding of data using sessions: '
                for scan_i in Sessions:
                    history_message += '{0}, '.format(scan_i)
        else:
            if not file_extension:
                file_extension='_sess{0}'.format(Sessions)
            history_message='Gridding of data using session {0}'.format(Sessions)

    if len(filelist) == 0:
        warnings.warn('There are no FITS files to process in '+rootdir+'/'+region+'/'+dirname)
        return
    # check that every file in the filelist is valid
    # If not then remove it and send warning message
    for file_i in filelist:
        try:
            fits.open(file_i)
        except:
            warnings.warn('file {0} is corrupted'.format(file_i))
            filelist.remove(file_i)
    #pull a test structure
    s = fits.getdata(filelist[0])

    c = 299792458.
    nu0 = s[0]['RESTFREQ']
    
    Data_Unit = s[0]['TUNIT7']
    beamSize = 1.22 * (c/nu0/100.0)*180/np.pi # in degrees

    naxis3 = len(s[0]['DATA'][startChannel:endChannel])
    nuslice = (range(naxis3))[startChannel:endChannel]

    # Default behavior is to park the object velocity at the center channel in the VRAD-LSR frame
    crval3 = s[0]['RESTFREQ']*(1-s[0]['VELOCITY']/c)
    crpix3 = s[0]['CRPIX1']-startChannel
    ctype3 = s[0]['CTYPE1']
    cdelt3 = s[0]['CDELT1']

    w = wcs.WCS(naxis=3)

    w.wcs.restfrq = nu0
    w.wcs.radesys = s[0]['RADESYS']
    w.wcs.equinox = s[0]['EQUINOX']
    w.wcs.specsys = 'LSRK'  # We are forcing this conversion to make nice cubes.
    w.wcs.ssysobs = 'TOPOCENT'
       
    if templateHeader is None:
        wcsdict = autoHeader(filelist, beamSize = beamSize, 
                             pixPerBeam = pixPerBeam)
        w.wcs.crpix = [wcsdict['CRPIX1'],wcsdict['CRPIX2'],crpix3]
        w.wcs.cdelt = np.array([wcsdict['CDELT1'],wcsdict['CDELT2'],cdelt3])
        w.wcs.crval = [wcsdict['CRVAL1'],wcsdict['CRVAL2'],crval3]
        w.wcs.ctype = [wcsdict['CTYPE1'],wcsdict['CTYPE2'],ctype3]
        naxis2 = wcsdict['NAXIS2']
        naxis1 = wcsdict['NAXIS1']
    else:
        w.wcs.crpix = [templateHeader['CRPIX1'],
                       templateHeader['CRPIX2'],crpix3]
        w.wcs.cdelt = np.array([templateHeader['CDELT1'],
                                templateHeader['CDELT2'],cdelt3])
        w.wcs.crval = [templateHeader['CRVAL1'],
                       templateHeader['CRVAL2'],crval3]
        w.wcs.ctype = [templateHeader['CTYPE1'],
                       templateHeader['CTYPE2'],ctype3]
        naxis2 = templateHeader['NAXIS2']
        naxis1 = templateHeader['NAXIS1']
    outCube = np.zeros((int(naxis3),int(naxis2),int(naxis1)))
    outWts = np.zeros((int(naxis2),int(naxis1)))

    xmat,ymat = np.meshgrid(np.arange(naxis1),np.arange(naxis2),indexing='ij')
    xmat = xmat.reshape(xmat.size)
    ymat = ymat.reshape(ymat.size)
    xmat = xmat.astype(np.int)
    ymat = ymat.astype(np.int)

    ctr = 0
    for thisfile in filelist:
        ctr+=1
        s = fits.open(thisfile)
        print("Now processing {0}".format(thisfile))
        print("This is file {0} of {1}".format(ctr,len(filelist)))

        nuindex = np.arange(len(s[1].data['DATA'][0]))

        EW = VelocityInterval(-5e3,5e3)
        EW.applyshift(8.500e3)

        for spectrum in console.ProgressBar((s[1].data)):
            # Generate Baseline regions
            if hasattr(baselineGenerator,'__call__'):
                baselineRegion = baselineGenerator(spectrum,EmissionWindow = EW)
                baselineIndex = np.concatenate([nuindex[ss] for ss in baselineRegion])
            else:
                baselineIndex = np.concatenate([nuindex[ss] for ss in baselineRegion])

            specData = spectrum['DATA']
            #baseline fit
            if doBaseline & np.all(np.isfinite(specData)):
                specData = baselineSpectrum(specData,order=blorder,
                                            baselineIndex=baselineIndex)

            # This part takes the TOPOCENTRIC frequency that is at
            # CRPIX1 (i.e., CRVAL1) and calculates the what frequency
            # that would have in the LSRK frame with freqShiftValue.
            # This then compares to the desired frequency CRVAL3.

            DeltaNu = freqShiftValue(spectrum['CRVAL1'],-spectrum['VFRAME'])-crval3
            DeltaChan = DeltaNu/cdelt3
            specData = channelShift(specData,-DeltaChan)
            outslice = (specData)[startChannel:endChannel]
            spectrum_wt = np.isfinite(outslice).astype(np.float)
            outslice = np.nan_to_num(outslice)
            xpoints,ypoints,zpoints = w.wcs_world2pix(spectrum['CRVAL2'],
                                                      spectrum['CRVAL3'],
                                                      spectrum['CRVAL1'],0)
            tsys = spectrum['TSYS']
            if (tsys>10) and (xpoints > 0) and (xpoints < naxis1) \
                    and (ypoints >0) and (ypoints < naxis2):
                pixelWeight,Index = gridFunction(xmat,ymat,
                                                 xpoints,ypoints,
                                                 pixPerBeam = pixPerBeam)
                vector = np.outer(outslice*spectrum_wt,pixelWeight/tsys**2)
                vector_wts = np.outer(spectrum_wt,pixelWeight/tsys**2)
                wts = pixelWeight/tsys**2
                outCube[:,ymat[Index],xmat[Index]] += vector
                outWts[ymat[Index],xmat[Index]] += wts
# Temporarily do a file write for every batch of scans.
        outWtsTemp = np.copy(outWts)
        outWtsTemp.shape = (1,)+outWtsTemp.shape
        outCubeTemp = np.copy(outCube)
        outCubeTemp /= outWtsTemp

        hdr = fits.Header(w.to_header())
        hdr = addHeader_nonStd( hdr, beamSize, Data_Unit)
        #
        hdu = fits.PrimaryHDU(outCubeTemp,header=hdr)
        hdu.writeto(dirname+file_extension+'.fits',clobber=True)

    outWts.shape = (1,)+outWts.shape
    outCube /= outWts

    # Create basic fits header from WCS structure
    hdr = fits.Header(w.to_header())
    # Add non standard fits keyword
    hdr = addHeader_nonStd( hdr, beamSize, Data_Unit)
    # Adds history message
    hdr.add_history(history_message)
    hdr.add_history('Using GAS pipeline version {0}'.format(__version__))    
    #
    hdu = fits.PrimaryHDU(outCube,header=hdr)
    hdu.writeto(dirname+file_extension+'.fits',clobber=True)

    w2 = w.dropaxis(2)
    hdr2 = fits.Header(w2.to_header())
    hdu2 = fits.PrimaryHDU(outWts,header=hdr2)
    hdu2.writeto(dirname+file_extension+'_wts.fits',clobber=True)

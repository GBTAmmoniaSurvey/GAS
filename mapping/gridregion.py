from sdpy import makecube
import numpy as np
import glob
from astropy.io import fits
import astropy.wcs as wcs
import itertools
from scipy.special import jv as besselJ
import pdb

import astropy.utils.console as console
rootdir = '/lustre/pipeline/scratch/GAS/NGC1333'
dirname = 'NGC1333_NH3_11'

filelist = glob.glob(rootdir+'/'+dirname+'/*feed1*sess1*fits')

def sincGrid(xpix,ypix,xdata,ydata, pixPerBeam = None):
    a = 1.55*(3.0/pixPerBeam)
    b = 2.52*(3.0/pixPerBeam)
    Rsup = 1.0
    dmin = 1e-4
    dx = (xdata - xpix)/pixPerBeam
    dy = (ydata - ypix)/pixPerBeam

    pia = np.pi/a
    b2 = 1./(b**2)
    distance = np.sqrt(dx**2+dy**2)
    
    ind  = (np.where(distance<=Rsup))
    d = distance[ind]
    wt = besselJ(d*pia,1)/\
        (pia*d)*\
        np.exp(-d**2*b2)
#    wt[ind] = np.exp(-distance[ind]**2*b2)*\
#              np.sin(pia*distance[ind])/\
#              (pia*distance[ind])
    wt[(d<dmin)]=1.0
    return(wt,ind)



#makecube.generate_header(restfreq='23.69447e9')

RAlist = []
DEClist = []

for thisfile in filelist:
    s = fits.getdata(thisfile)
    RAlist = RAlist + [s['CRVAL2']]
    DEClist = DEClist +[s['CRVAL3']]

c = 299792458.
nu0 = s[0]['RESTFREQ']
beamSize = 1.22 * (c/nu0/100.0)*180/np.pi # in arcsec
pixPerBeam = 3.0

startChan = 1024
endChan = 4096/4*3

longitude = np.array(list(itertools.chain(*RAlist)))
latitude = np.array(list(itertools.chain(*DEClist)))

minLon = longitude.min()
maxLon = longitude.max()
minLat = latitude.min()
maxLat = latitude.max()


naxis3 = len(s[0]['DATA'])
nuslice = (range(naxis3))[startChan:endChan]
naxis3 = len(nuslice)

crpix3 = s[0]['CRPIX1']-startChan
crval3 = s[0]['CRVAL1']
ctype3 = s[0]['CTYPE1']
cdelt3 = s[0]['CDELT1']

naxis2 = np.ceil((maxLat-minLat)/(beamSize/pixPerBeam)+2*pixPerBeam)
crpix2 = naxis2/2
cdelt2 = beamSize/pixPerBeam
crval2 = (maxLat+minLat)/2
ctype2 = s[0]['CTYPE3']+'--TAN'
cdelt1 = beamSize/pixPerBeam
naxis1 = np.ceil((maxLon-minLon)/(beamSize/pixPerBeam)*\
                     np.cos(crval2/180*np.pi)+2*pixPerBeam)
crpix1 = naxis1/2
crval1 = (minLon+maxLon)/2
ctype1 = s[0]['CTYPE2']+'---TAN'

w = wcs.WCS(naxis=3)
w.wcs.crpix = [crpix1,crpix2,crpix3]
w.wcs.cdelt = np.array([cdelt1,cdelt2,cdelt3])
w.wcs.crval = [crval1,crval2,crval3]
w.wcs.ctype = [ctype1,ctype2,ctype3]

outCube = np.zeros((naxis3,naxis2,naxis1))
outWts = np.zeros((naxis2,naxis1))

xmat,ymat = np.meshgrid(np.arange(naxis1),np.arange(naxis2),indexing='ij')
xmat = xmat.reshape(xmat.size)
ymat = ymat.reshape(ymat.size)
xmat = xmat.astype(np.int)
ymat = ymat.astype(np.int)

gridFunction = sincGrid

ctr = 0
for thisfile in filelist:
    ctr+=1
    s = fits.open(thisfile)
    print("Now processing {0}".format(thisfile))
    print("This is file {0} of {1}".format(ctr,len(filelist)))
    for spectrum in console.ProgressBar((s[1].data)):
        outslice = (spectrum['DATA'])[nuslice]
        spectrum_wt = np.isfinite(outslice).astype(np.float)
        outslice = np.nan_to_num(outslice)
        xpoints,ypoints,zpoints = w.wcs_world2pix(spectrum['CRVAL2'],
                                                  spectrum['CRVAL3'],
                                                  spectrum['CRVAL1'],0)
        pixelWeight,Index = gridFunction(xmat,ymat,
                                            xpoints,ypoints,
                                            pixPerBeam = 3.0)
        vector = np.outer(outslice*spectrum_wt,pixelWeight/spectrum['TSYS'])
        vector_wts = np.outer(spectrum_wt,pixelWeight/spectrum['TSYS'])
        wts = pixelWeight/spectrum['TSYS']
        outCube[:,ymat[Index],xmat[Index]] += vector
        outWts[ymat[Index],xmat[Index]] += wts

outWts.shape = (1,)+outWts.shape
outCube /= outWts

hdr = fits.Header(w.to_header())
hdu = fits.PrimaryHDU(outCube,header=hdr)
hdu.writeto(dirname+'.fits',clobber=True)

w2 = w.dropaxis(2)
hdr2 = fits.Header(w2.to_header())
hdu2 = fits.PrimaryHDU(outWts,header=hdr2)
hdu2.writeto(dirname+'_wts.fits',clobber=True)

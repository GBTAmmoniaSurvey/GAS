
import astropy
import numpy as np
import scipy
import astropy.wcs as wcs
import pdb
import astropy.io.fits as fits

def sincGrid(xpix,ypix,xdata,ydata, pixPerBeam = None):
    a = 1.55
    b = 2.52
    Rsup = 3.
    dmin = 1e-4
    dx = (xdata - xpix)/pixPerBeam
    dy = (ydata - ypix)/pixPerBeam

    pia = np.pi/a
    b2 = 1./(b**2)
    distance = np.sqrt(dx**2+dy**2)
    wt = np.zeros(len(distance))
    
    ind  = (np.where(distance<=Rsup))
    wt[ind] = np.exp(-distance[ind]**2*b2)*\
              np.sin(pia*distance[ind])/\
              (pia*distance[ind])
    wt[(distance<dmin)]=1.0
    return(wt,ind)

def griddata(inFile, beamSize = 0.0087, pixPerBeam = 3.0,
             buffer = 1.1, wcsObject = None,
             gridFunction = sincGrid):

    hdu = fits.open(inFile)
    data = hdu[1].data
    longitude = data['CRVAL2']
    latitude = data['CRVAL3']
    minLon = longitude.min()
    maxLon = longitude.max()
    minLat = latitude.min()
    maxLat = latitude.max()

    naxis3 = len(data[0]['DATA'])
    crpix3 = data[0]['CRPIX1']
    crval3 = data[0]['CRVAL1']
    ctype3 = data[0]['CTYPE1']
    cdelt3 = data[0]['CDELT1']

    naxis2 = np.ceil((maxLat-minLat)/(beamSize/pixPerBeam)+2*pixPerBeam)
    crpix2 = naxis2/2
    cdelt2 = beamSize/pixPerBeam
    crval2 = (maxLat+minLat)/2
    ctype2 = data[0]['CTYPE3']+'--TAN'
# The projections need better definition

    cdelt1 = beamSize/pixPerBeam
    naxis1 = np.ceil((maxLon-minLon)/(beamSize/pixPerBeam)*\
                     np.cos(crval2/180*np.pi)+2*pixPerBeam)
    crpix1 = naxis1/2
    crval1 = (minLon+maxLon)/2
    ctype1 = data[0]['CTYPE2']+'---TAN'
    w = wcs.WCS(naxis=3)
    w.wcs.crpix = [crpix1,crpix2,crpix3]
    w.wcs.cdelt = np.array([cdelt1,cdelt2,cdelt3])
    w.wcs.crval = [crval1,crval2,crval3]
    w.wcs.ctype = [ctype1,ctype2,ctype3]

    xpoints,ypoints,zpoints = w.wcs_world2pix(data['CRVAL2'],
                                              data['CRVAL3'],
                                              data['CRVAL1'],0)
    dataPlane = data['DATA']
    outCube = np.zeros((naxis1,naxis2,naxis3))
    for ii in np.arange(naxis1):
        for jj in np.arange(naxis2):
            pixelWeight,useIndex = gridFunction(ii,jj,
                                                xpoints,ypoints,
                                                pixPerBeam = pixPerBeam)
 #           print(ii,jj)
            if useIndex and pixelWeight.sum()>0.3:
                spectrum = np.zeros(naxis3)
                for kk,idx in enumerate(useIndex[0]):
                    if np.all(np.isfinite(dataPlane[idx,:])):
                        spectrum = spectrum+pixelWeight[idx]*\
                                   dataPlane[idx,:]
                outCube[ii,jj,:]=spectrum/pixelWeight.sum()
    return(outCube)
            
        



    
    
       

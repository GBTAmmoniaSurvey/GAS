import astropy.coordinates as coords
import catalogs

import numpy as np


def VlsrByCoord(RA, Dec, region = 'OrionA', regionCatalog = None):

    if 'NGC1333' in region:
        v0 = 7.6687698945600005\
            -0.10015533*(RA-52.2275368656)\
            +1.29725507*(Dec-31.2380245486)
        return(v0)

    if 'OrionA' in region:
        coeffs = [-2.8256074 , -4.65791997,  9.14502305]
        v0 = coeffs[2]+\
            coeffs[0]*(RA-83.446122802665869)+\
            coeffs[1]*(Dec+6.0050560818354661)
        return(v0)

    if regionCatalog is None:
        regionCatalog = catalogs.GenerateRegions()

    try:
        v0 = regionCatalog[regionCatalog['Region name']==
                           region]['VLSR'].data.data[0]
        return v0
    except IndexError:
        return np.nan

def FitGradient(vlsrData, vlsrWCS,region='OrionA'):
    
    ygood,xgood = np.where((vlsrData!=0)*np.isfinite(vlsrData))
    ra, dec = vlsrWCS.wcs_pix2world(xgood,ygood,0)
    rabar= np.median(ra)
    decbar = np.median(dec)
    vbar = np.median(vlsrData[ygood,xgood])
    vvals = vlsrData[ygood,xgood]-vbar

    A = np.r_[[np.ones(len(ra))],[ra-rabar],[dec-decbar]]
    invcov = np.linalg.inv(np.matrix(A)*np.matrix(A.T))
    
    w = np.matrix([np.sum(vvals),
                   np.sum(vvals*(ra-rabar)),
                   np.sum(vvals*(dec-decbar))]).T
    coeffs = np.array(invcov*w)

    stringarr = ['Gradient Values',
                 'meanRA  = {0}'.format(rabar),
                 'meanDec = {0}'.format(decbar),
                 'gradRA  = {0}'.format(coeffs[1]),
                 'gradDec = {0}'.format(coeffs[2]),
                 'voffset = {0}'.format(coeffs[0]+vbar)]
    for i in stringarr:
        print(i)

    return

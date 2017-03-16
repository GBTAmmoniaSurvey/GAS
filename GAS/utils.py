import catalogs
import numpy as np


def VlsrByCoord(RA, Dec, region='OrionA', regionCatalog=None):

    """
    For the GAS regions determine the appropriate v0 given coordinates.
    If a velocity gradient has been fit, it is hard-coded here.  Otherwise,
    this uses the catalog velocities.

    Parameters
    ----------
    RA : float
        The right ascension of the spectrum
    Dec: float
        The declination of the spectrum
    Region : str
        The GAS name of the region, used in catalog lookups.
    regionCatalog : astropy.table.Table
        Table containing the catalog values of the velocity to prevent
        repeated catalog parsing

    Returns
    -------
    v0 : float
        Estimate of central velocity of emission spectrum in km/s
    """

    if 'NGC1333' in region:
        v0 = 7.6687698945600005\
            - 0.10015533 * (RA - 52.2275368656)\
            + 1.29725507 * (Dec - 31.2380245486)
        return(v0)

    if 'OrionA' in region:
        coeffs = [3.1792362, 2.36881382, 9.25135413]
        v0 = coeffs[2] + \
            coeffs[0] * (RA - 83.8214130221) + \
            coeffs[1] * (Dec + 5.60497968304)
        return(v0)

    if regionCatalog is None:
        regionCatalog = catalogs.GenerateRegions()

    try:
        v0 = regionCatalog[regionCatalog['Region name'] ==
                           region]['VLSR'].data.data[0]
        return v0
    except IndexError:
        return np.nan


def FitGradient(vlsrData, vlsrWCS):

    """
    Utility function for fitting a velocity gradient.  Gives coefficients to
    hardwire into the VlsrByCoord function.

    Parameters
    ----------
    vlsrData : np.array
        Image of the vlsr data
    vlsrWCS : astropy.wcs.WCS
        WCS object corresponding to the image

    Returns
    -------
    None -- output print to screen.
    """

    ygood, xgood = np.where((vlsrData != 0) * np.isfinite(vlsrData))
    ra, dec = vlsrWCS.wcs_pix2world(xgood, ygood, 0)
    rabar = np.median(ra)
    decbar = np.median(dec)
    vbar = np.median(vlsrData[ygood, xgood])
    vvals = vlsrData[ygood, xgood] - vbar

    A = np.r_[[np.ones(len(ra))], [ra - rabar], [dec - decbar]]
    invcov = np.linalg.inv(np.matrix(A) * np.matrix(A.T))

    w = np.matrix([np.sum(vvals),
                   np.sum(vvals * (ra - rabar)),
                   np.sum(vvals * (dec - decbar))]).T
    coeffs = np.array(invcov * w)

    stringarr = ['Gradient Values',
                 'meanRA  = {0}'.format(rabar),
                 'meanDec = {0}'.format(decbar),
                 'gradRA  = {0}'.format(coeffs[1]),
                 'gradDec = {0}'.format(coeffs[2]),
                 'voffset = {0}'.format(coeffs[0] + vbar)]
    for i in stringarr:
        print(i)

    return

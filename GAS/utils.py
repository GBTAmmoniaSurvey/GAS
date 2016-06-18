import astropy.coordinates as coords
import catalogs

import numpy as np


def VlsrByCoord(RA, Dec, region = 'OrionA', regionCatalog = None):

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

def FitGradient(vlsr_file,region='OrionA'):
    raise NotImplementedError("Coming soon to a codebase near you")
    return

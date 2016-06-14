import astropy.coordinates as coords
import catalogs


def VlsrByCoord(RA, Dec, region = 'OrionA'):

    regions = catalogs.GenerateRegions()

    if 'OrionA' in region:
        coeffs = [-2.8256074 , -4.65791997,  9.14502305]
        v0 = coeffs[2]+\
            coeffs[0]*(RA-83.446122802665869)+\
            coeffs[1]*(Dec+6.0050560818354661)

    return(v0)

def FitGradient(vlsr_file,region='OrionA'):
    
    return coeffs, ra0, dec0

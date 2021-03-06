import pyspeckit
import astropy.io.fits as fits
import numpy as np
import os
from .first_look import trim_edge_cube
from spectral_cube import SpectralCube
from pyspeckit.spectrum.models.ammonia_constants import voff_lines_dict
import astropy.constants as con
import astropy.units as u
import matplotlib.pyplot as plt
import aplpy
from skimage.morphology import remove_small_objects,closing,disk,opening
from .gauss_fit import gauss_fitter
from .run_first_look import trim_edge_spectral_cube
import glob
from . import catalogs

from pyspeckit.spectrum.models import ammonia
from .config import plottingDictionary

'''
Set of scripts to update and streamline masks based on updated DR2 
moment and property maps
'''


def mask_all_data(regions=None,file_extension='all_rebase3',release='all'):
    if regions is None:
        RegionCatalog = catalogs.GenerateRegions(release=release)
    else:
        RegionCatalog = catalogs.GenerateRegions(release=release)
        keep = [idx for idx, row in enumerate(RegionCatalog) if row['Region name'] in regions]
        RegionCatalog = RegionCatalog[keep]
        
    for ThisRegion in RegionCatalog:
        region=ThisRegion['Region name']
        trim_cubes(region=region,file_extension=file_extension)
        flag_all_data(region=region,file_extension=file_extension)

def mask_binary(imageHDU,LowestContour=3,selem=np.array([[0,1,0],[1,1,1],[0,1,0]])):
    """
    Way to mask 'island pixels' in property maps
    """
    from scipy.ndimage import binary_opening
    imageMap = imageHDU.data
    mask = binary_opening(imageMap > LowestContour, selem)
    #imageHDU.close()
    #MaskedMap = mask*imageMap
    #imageHDU[0].data = MaskedMap
    return mask

def trim_cubes(region='OrionA',file_extension='all_rebase3',propertyMaps=True):
    '''
    Trim cube edges beyond initial trimming in earlier pipeline steps. 
    Since trim looks for edges in the data, CANNOT use on already flagged data.
    Use updated moment map as mask for other files
    '''
    root = file_extension

    # Cubes:
    line_list=['NH3_11','NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    #line_list = ['NH3_11','NH3_22']
    for line in line_list:
        # Moment first to create mask
        try:
            moment = fits.open('{0}/{0}_{1}_{2}_mom0_QA.fits'.format(region,line,file_extension))
            moment_data = moment[0].data
            moment_hdr  = moment[0].header
            moment.close()
            trim_edge_cube(moment_data)
            mask = np.isfinite(moment_data)
            # Write out new moment
            fits.writeto('{0}/{0}_{1}_{2}_mom0_QA_trim.fits'.format(region,line,file_extension),
                         moment_data,moment_hdr,overwrite=True)
            # Next cubes
            filein = '{0}/{0}_{1}_{2}.fits'.format(region,line,file_extension)
            # trim_edge_cube doesn't work on a spectral cube object. 
            cube = SpectralCube.read(filein)
            cube2 = cube.with_mask(mask)
            cube2.write('{0}/{0}_{1}_{2}_trim.fits'.format(region,line,file_extension),overwrite=True)
            # And rms. Note that QA rms file has interior masked regions for Orion A, others?
            rms = fits.open('{0}/{0}_{1}_{2}_rms_QA.fits'.format(region,line,file_extension))
            rms_data = rms[0].data
            rms_hdr  = rms[0].header
            rms.close()
            trim_edge_cube(rms_data)
            fits.writeto('{0}/{0}_{1}_{2}_rms_QA_trim.fits'.format(region,line,file_extension),
                         rms_data,rms_hdr,overwrite=True)
        except:
            print '{0}/{0}_{1}_{2}_mom0_QA.fits does not exist.'.format(region,line,file_extension)
    # Use NH3 (1,1) moment map as mask for property map
    if propertyMaps:
        fit_file = '{0}/{0}_parameter_maps_{1}.fits'.format(region,file_extension)
        moment = fits.open('{0}/{0}_NH3_11_{1}_mom0_QA_trim.fits'.format(region,file_extension))
        moment_data = moment[0].data
        propMap = fits.open(fit_file)
        propMap_data = propMap[0].data
        propMap_hdr  = propMap[0].header
        propMap.close()    
        for plane_i in range(len(propMap_data)):
            data_i = propMap_data[plane_i,:,:]
            data_i[~np.isfinite(moment_data)] = np.nan
            propMap_data[plane_i,:,:] = data_i
        fits.writeto('{0}/{0}_parameter_maps_{1}_trim.fits'.format(region,file_extension),
                     propMap_data,propMap_hdr,overwrite=True)
    

def flag_all_data(region='OrionA',file_extension='all_rebase3'):
    """
    Flag cubefit results based on S/N in integrated intensity.
    Also flag poorly constrained fits (where Tk, Tex hit minimum values)
    Outputs individual .fits files for parameters and uncertainties: Tkin, Tex, Vc, sigmaV, NNH3
    Remove moment 0 flagging as update_NH3_moment + edge masking is better
    Update (16/9/13): Include masking of 'island pixels' based on S/N in NH3 (1,1) moment map
                      Set flagged values to zero rather than NaN
                      Output masked cube file
                      Cleaned up
    Update (19/7/12): Removed S/N limit on (2,2) line, base masking of, e.g., Tk only on 
                      uncertainties in fitted parameters
    Update (19/8/9):  Additional mask on N, Tk where no Tex values, and upper limit on e(N)

    Parameters
    ----------
    region : str
        Name of region to reduce
    blorder : int
        order of baseline removed
    file_extension: : str
        filename
    """
    import matplotlib.pyplot as plt
    import aplpy

    flagMinTex = 2.8
    flagMaxeTex = 3.0

    flagMinTk = 5.0
    flagMaxeTk = 3.0

    flagMaxeN = 0.15

    flagSN11 = 3.0

    root = file_extension
    if not os.path.exists('{0}/parameterMaps'.format(region)):
        os.mkdir('{0}/parameterMaps'.format(region))
        
    hdu=fits.open("{0}/{0}_parameter_maps_{1}_trim.fits".format(region,root))
    hd_cube=hdu[0].header
    cube=hdu[0].data
    hdu.close()

    rms11hdu = fits.open("{0}/{0}_NH3_11_{1}_rms_QA.fits".format(region,root))
    rms11data = rms11hdu[0].data
    rms11hdu.close()
    m0_11 = fits.open("{0}/{0}_NH3_11_{1}_mom0_QA.fits".format(region,root))
    m0_11data = m0_11[0].data
    hd11 = m0_11[0].header
    m0_11.close()

    sn11 = m0_11data/rms11data
    sn11HDU = fits.PrimaryHDU(sn11,hd11)

    # Get binary mask
    selem=np.array([[0,1,0],[1,1,1],[0,1,0]])
    LowestContour=3
    pixel_mask = mask_binary(sn11HDU,LowestContour=LowestContour,selem=selem)

    # Get Tex, Tk files for mask
    tex  = np.copy(cube[1,:,:])
    etex = np.copy(cube[7,:,:])
    tk   = np.copy(cube[0,:,:])
    etk  = np.copy(cube[6,:,:])
    eN   = np.copy(cube[8,:,:])

    hd = hd_cube.copy()
    rm_key=['NAXIS3','CRPIX3','CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3']
    for key_i in rm_key:
        hd.remove(key_i)
    hd['NAXIS']= 2
    hd['WCSAXES']= 2

    # Tkin
    hd['BUNIT']='K'
    param=cube[0,:,:]
    eparam = cube[6,:,:]
    parMask = ((tex >flagMinTex) & (param >flagMinTk) & \
               (etex < flagMaxeTex) & (eparam < flagMaxeTk) & (eparam != 0) & (eN<flagMaxeN) & (eN > 0))
    param = param * pixel_mask * parMask
    eparam = eparam * pixel_mask * parMask
    file_out="{0}/parameterMaps/{0}_NH3_Tkin_{1}_masked.fits".format(region,root)
    fits.writeto(file_out, param, hd, overwrite=True)
    # Write out uncertainties
    file_out="{0}/parameterMaps/{0}_NH3_eTkin_{1}_masked.fits".format(region,root)
    fits.writeto(file_out, eparam, hd, overwrite=True)
    # Update cube
    cube[0,:,:] = param
    cube[6,:,:] = eparam

    #Tex
    hd['BUNIT']='K'
    param=cube[1,:,:]
    eparam=cube[7,:,:]
    parMask = ((sn11>=flagSN11) & (tex > flagMinTex) & (tk > flagMinTk) & (etex < flagMaxeTex))
    param = param * pixel_mask * parMask
    eparam = eparam * pixel_mask * parMask
    file_out="{0}/parameterMaps/{0}_NH3_Tex_{1}_masked.fits".format(region,root)
    fits.writeto(file_out, param, hd, overwrite=True)
    # Write out uncertainties
    file_out="{0}/parameterMaps/{0}_NH3_eTex_{1}_masked.fits".format(region,root)
    fits.writeto(file_out, eparam, hd, overwrite=True)
    # Update cube
    cube[1,:,:] = param
    cube[7,:,:] = eparam

    # N_NH3
    hd['BUNIT']='cm-2'
    param=cube[2,:,:]
    eparam=cube[8,:,:]
    parMask = ((sn11>=flagSN11) & (tex > flagMinTex) & (tk > flagMinTk) & (etex < flagMaxeTex) &
               (eparam < flagMaxeN) & (eparam > 0))
    param = param * pixel_mask * parMask
    eparam = eparam * pixel_mask * parMask
    file_out="{0}/parameterMaps/{0}_NH3_N_{1}_masked.fits".format(region,root)
    fits.writeto(file_out, param, hd, overwrite=True)
    # Write out uncertainties
    file_out="{0}/parameterMaps/{0}_NH3_eN_{1}_masked.fits".format(region,root)
    fits.writeto(file_out, eparam, hd, overwrite=True)    # sigma
    # Update cube
    cube[2,:,:] = param
    cube[8,:,:] = eparam

    # Use same flags for Vlsr, sigma
    # Sigma
    hd['BUNIT']='km/s'
    param=cube[3,:,:]
    eparam=cube[9,:,:]
    parMask = ((sn11 >= flagSN11))
    param = param * pixel_mask * parMask
    eparam = eparam * pixel_mask * parMask
    file_out="{0}/parameterMaps/{0}_NH3_Sigma_{1}_masked.fits".format(region,root)
    fits.writeto(file_out, param, hd, overwrite=True)
    # Write out uncertainties
    file_out="{0}/parameterMaps/{0}_NH3_eSigma_{1}_masked.fits".format(region,root)
    fits.writeto(file_out, eparam, hd, overwrite=True)
    # Update cube
    cube[3,:,:] = param
    cube[9,:,:] = eparam

    # Vlsr
    hd['BUNIT']='km/s'
    param=cube[4,:,:]
    eparam=cube[10,:,:]
    parMask = ((sn11 >= flagSN11))
    param = param * pixel_mask * parMask
    eparam = eparam * pixel_mask * parMask
    file_out="{0}/parameterMaps/{0}_NH3_Vlsr_{1}_masked.fits".format(region,root)
    fits.writeto(file_out, param, hd, overwrite=True)
    # Write out uncertainties
    file_out="{0}/parameterMaps/{0}_NH3_eVlsr_{1}_masked.fits".format(region,root)
    fits.writeto(file_out, eparam, hd, overwrite=True)
    # Update cube
    cube[4,:,:] = param
    cube[10,:,:] = eparam

    # Write out new cube
    cube_out="{0}/{0}_parameter_maps_{1}_masked.fits".format(region,root)
    fits.writeto(cube_out,cube,hd_cube,overwrite=True)

import pyspeckit
import astropy.io.fits as fits
import numpy as np
import os
import warnings
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
from . import catalogs
import glob

from pyspeckit.spectrum.models import ammonia
from .config import plottingDictionary

def update_NH3_moment0(region_name='L1688', file_extension='DR1_rebase3', threshold=0.0125, trim_edge=True, save_masked=False, save_model=False):
    """
    Function to update moment calculation based on centroid velocity from line fit.
    For a given NH3(1,1) cube, we check which channels have flux in the model cube, 
    and then use those channels as the appropiate channels for integration.

    Based on code provided by Vlas Sokolov

    Parameters
    ----------
    region : str
        Name of region to re-calculate moment map
    file_extension : str
        filename extension
    threshold : float
        minimum threshold in model cube used to identify channels with emission
        Default is down to the machine precision, a better result could be 
        obtained with 0.0125
    trim_edge : Boolean
        Keyword to trim noisy edges of output maps using disk erode 
    save_masked : Boolean
        Keyword to store the masked cube used in the integrated intensity calculation.
        This is useful to 

    Usage: 
    import GAS
    GAS.PropertyMaps.update_NH3_moment0(region_name='NGC1333', file_extension='DR1_rebase3', threshold=0.0125, save_masked=True)

    """
    fit_file='{0}/{0}_parameter_maps_{1}.fits'.format(region_name,file_extension)
    for line_i in ['11','22']:
        file_in ='{0}/{0}_NH3_{2}_{1}.fits'.format(region_name,file_extension,line_i)
        file_out='{0}/{0}_NH3_{2}_{1}_mom0_QA.fits'.format(region_name,file_extension,line_i)
        file_rms_init = '{0}/{0}_NH3_{2}_{1}_rms.fits'.format(region_name,file_extension,line_i)
        file_rms='{0}/{0}_NH3_{2}_{1}_rms_QA.fits'.format(region_name,file_extension,line_i)
        file_rms_mom='{0}/{0}_NH3_{2}_{1}_mom0_sigma_QA.fits'.format(region_name,file_extension,line_i)
        file_temp='{0}/{0}_NH3_{2}_{1}_masked_temp.fits'.format(region_name,file_extension,line_i)
        plot_param = plottingDictionary[region_name]
        # Load pyspeckit cube
        pycube = pyspeckit.Cube(file_in)
        rms_init = fits.getdata(file_rms_init)
        if 'FITTYPE' in fits.getheader(fit_file):
            # 'FITTYPE' is not present in old versions of the parameter files
            pycube.load_model_fit( fit_file, npars=6, npeaks=1)
        else:
    	    if not 'cold_ammonia' in pycube.specfit.Registry.multifitters:
                pycube.specfit.Registry.add_fitter('cold_ammonia',ammonia.cold_ammonia_model(),6)
            pycube.load_model_fit( fit_file, npars=6, npeaks=1, fittype='cold_ammonia')
        # If threshold is not defined, then use the machine accuracy
        if threshold == None:
            threshold=np.finfo(pycube.data.dtype).eps
        # Get model cube from pyspeckit, this take some time
        modelcube = pycube.get_modelcube()
        # Use spectral cube to calculate integrated intensity maps
        cube_raw = SpectralCube.read(file_in)
        # in km/s not Hz
        cube = cube_raw.with_spectral_unit(u.km / u.s,velocity_convention='radio')
        if trim_edge:
            cube = trim_edge_spectral_cube(cube)
        vaxis=cube.spectral_axis
        dv=np.abs(vaxis[1]-vaxis[0])
        if save_model:
            model_scube = SpectralCube(data=modelcube,wcs=cube.wcs)
            model_scube.write('{0}/{0}_NH3_{2}_{1}_model.fits'.format(region_name,file_extension,line_i),overwrite=True)
        # define mask 
        mask3d = modelcube > threshold
        # What to do with pixels without signal
        # Calculate mean velocity and velocity dispersion
        vmap=pycube.parcube[4,:,:]
        sigma_map=pycube.parcube[3,:,:]
        vmean=np.mean(vmap[vmap != 0])*u.km/u.s
        if line_i == '11':
            voff11 = voff_lines_dict['oneone']
            sigma_v=( np.median(sigma_map[vmap != 0]) + 0.15)*u.km/u.s
            total_spc = np.ones(cube.shape[0],dtype=np.bool)
            for deltav in voff11:
                total_spc *= (np.abs(cube.spectral_axis-(deltav*u.km/u.s+vmean))) > plot_param['sigma_mult'] * sigma_v
            total_spc = ~total_spc
        else:
            sigma_v=( np.median(sigma_map[vmap != 0]))*u.km/u.s
            total_spc=np.sqrt( (vaxis-vmean)**2)/sigma_v < plot_param['sigma_mult']
        # 
        im_mask=np.sum(mask3d, axis=0)
        # Here checking for bad fits:
        # No fits. Places where parameter uncertainties are zero or unreasonably low
        # Where Tk, Tex hit minimum allowed values
        # Fitted line width < 3 * error in line width
        # Added check on model fit amplitude vs. spectrum rms
        for ii in np.arange( im_mask.shape[1]):
            for jj in np.arange( im_mask.shape[0]):
                if ((im_mask[jj,ii] == 0) or (pycube.parcube[3,jj,ii] < 3*pycube.errcube[3,jj,ii]) or
                    (pycube.errcube[4,jj,ii] == 0) or (pycube.errcube[0,jj,ii] > 5) or
                    (pycube.errcube[1,jj,ii] < 0.01) or 
                    (max(modelcube[:,jj,ii]) < (3.*rms_init[jj,ii])) or 
		    (pycube.parcube[0,jj,ii] == 5) or (pycube.parcube[1,jj,ii] < 3)):
                    # Find vlsr of nearby fits
                    vmap_loc = vmap[jj-20:jj+20,ii-20:ii+20]
                    vmean_loc = np.mean(vmap_loc[vmap_loc != 0])*u.km/u.s
                    if np.isfinite(vmean_loc):
                        vshift = np.int((vmean - vmean_loc)/dv)
                        # Shift total_spc to match local vlsr
                        total_spc_shift = np.roll(total_spc,vshift)
                    else:
                        # Try bigger window
                        vmap_loc = vmap[jj-50:jj+50,ii-50:ii+50]
                        vmean_loc = np.mean(vmap_loc[vmap_loc != 0])*u.km/u.s
                        if np.isfinite(vmean_loc):
                            vshift = np.int((vmean - vmean_loc)/dv)
                            # Shift total_spc to match local vlsr
                            total_spc_shift = np.roll(total_spc,vshift)
                        else:
                            # Use original window if no good fits nearby
                            total_spc_shift = total_spc
                    mask3d[:,jj,ii] = total_spc_shift
        n_chan=np.sum(mask3d, axis=0)
        # create masked cube
        cube2 = cube.with_mask(mask3d)
        cube3 = cube.with_mask(~mask3d)
        #
        if save_masked:
            cube2.write( file_temp, overwrite=True)
        # calculate moment map
        moment_0 = cube2.moment(axis=0)
        moment_0.write( file_out, overwrite=True)
        rms=cube3.std(axis=0)
        rms.write( file_rms, overwrite=True)
        mom_0_rms=rms * dv * np.sqrt(n_chan)
        mom_0_rms.write( file_rms_mom, overwrite=True)

def update_rest_moment0(region_name='L1688', file_extension='DR1_rebase3', v_mean=2.5, sigma_v=0.3):
    """
    Update moment calculation for non-NH3 lines. Don't base on NH3 fit; these lines may be
    offset in velocity. Right now including velocity ranges by hand. 
    """
    for line in ['HC5N','HC7N_21_20','HC7N_22_21','C2S','NH3_33']:
        file_in ='{0}/{0}_{2}_{1}.fits'.format(region_name,file_extension,line)
        file_out='{0}/{0}_{2}_{1}_mom0_QA.fits'.format(region_name,file_extension,line)
        file_rms='{0}/{0}_{2}_{1}_rms_QA.fits'.format(region_name,file_extension,line)
        # Load pyspeckit cube
        pycube = pyspeckit.Cube(file_in)
        # Use spectral cube to calculate integrated intensity maps
        cube_raw = SpectralCube.read(file_in)
        # in km/s not Hz
        cube = cube_raw.with_spectral_unit(u.km / u.s,velocity_convention='radio')
        vaxis=cube.spectral_axis
        dv=np.abs(vaxis[1]-vaxis[0])
        vmean = v_mean*u.km/u.s
        sigmav = sigma_v*u.km/u.s
        total_spc=np.sqrt( (vaxis-vmean)**2)/sigmav < 3
        mask3d = np.ones(cube.shape,dtype=bool)
        for ii in np.arange(cube.shape[2]):
            for jj in np.arange(cube.shape[1]):
                mask3d[:,jj,ii] = total_spc
        # create masked cube
        cube2 = cube.with_mask(mask3d)
        cube3 = cube.with_mask(~mask3d)
        # calculate moment map
        moment_0 = cube2.moment(axis=0)
        moment_0.write( file_out, overwrite=True)
        rms=cube3.std(axis=0)
        rms.write( file_rms, overwrite=True)

def update_rest_moment0_2(region_name='L1688', file_extension='all_rebase3',
                          threshold=0.0125, save_masked=False,trim_edge=True,
                          vsys=5.*u.km/u.s,throw=2.*u.km/u.s):
    """
    Function to update moment calculation based on centroid velocity from line fit.
    For a given line cube, we check which channels have flux in the model cube, 
    and then use those channels as the appropiate channels for integration.

    Based on code provided by Vlas Sokolov

    Parameters
    ----------
    region : str
        Name of region to re-calculate moment map
    file_extension : str
        filename extension
    threshold : float
        minimum threshold in model cube used to identify channels with emission
        Default is down to the machine precision, a better result could be 
        obtained with 0.0125
    save_masked : Boolean
        Keyword to store the masked cube used in the integrated intensity calculation.
        This is useful to 
    trim_edge : Boolean
        Trim map edges (3 pixels in)
    vsys, throw: float w. units
        Describes window in velocity to be used for moment and rms maps 
        if no good fits to lines are found in the cube. 
    Usage: 
    import GAS
    GAS.PropertyMaps.update_rest_moment0_2(region_name='NGC1333', file_extension='DR1_rebase3', threshold=0.0125, save_masked=True)

    """
    
    for line in ['HC5N','HC7N_21_20','HC7N_22_21','C2S','NH3_33']:
        fit_file='{0}/{0}_{1}_{2}_param_cube_masked.fits'.format(region_name,line,file_extension)
        fit_model_file='{0}/{0}_{1}_{2}_gauss_cube.fits'.format(region_name,line,file_extension)
        file_in ='{0}/{0}_{1}_{2}.fits'.format(region_name,line,file_extension)
        file_out='{0}/{0}_{1}_{2}_mom0_QA.fits'.format(region_name,line,file_extension)
        file_rms='{0}/{0}_{1}_{2}_rms_QA.fits'.format(region_name,line,file_extension)
        file_rms_mom='{0}/{0}_{1}_{2}_mom0_sigma_QA.fits'.format(region_name,line,file_extension)
        file_temp='{0}/{0}_{1}_{2}_masked_temp.fits'.format(region_name,line,file_extension)
        # Check that files exist
        if os.path.isfile(file_in):
            # Load pyspeckit cube
            pycube = pyspeckit.Cube(file_in)
            pycube.load_model_fit( fit_file, npars=3, npeaks=1, fittype='gaussian')
            # If threshold is not defined, then use the machine accuracy
            if threshold == None:
                threshold=np.finfo(pycube.data.dtype).eps
            # Get model cube from pyspeckit. Is completely zero for Gaussian fits. Why?
            try:
            #modelcube = pycube.get_modelcube()
            # Using output model file from the Gaussian fitter instead
                model = pyspeckit.Cube(fit_model_file)
                modelcube = model.cube
                # Use spectral cube to calculate integrated intensity maps
                cube_raw = SpectralCube.read(file_in)
                # in km/s not Hz
                cube = cube_raw.with_spectral_unit(u.km / u.s,velocity_convention='radio')
                if trim_edge:
                    cube = trim_edge_spectral_cube(cube) 
                vaxis=cube.spectral_axis
                dv=np.abs(vaxis[1]-vaxis[0])
                # define mask 
                mask3d = modelcube > threshold
                # What to do with pixels without signal
                # Calculate mean velocity and velocity dispersion
                vmap=pycube.parcube[1,:,:]
                sigma_map=pycube.parcube[2,:,:]
                vmean=np.nanmean(vmap[vmap != 0])*u.km/u.s
                sigma_v=( np.nanmedian(sigma_map[vmap != 0]))*u.km/u.s
                if np.isfinite(vmean):
                    total_spc=np.sqrt( (vaxis-vmean)**2)/sigma_v < 3.0
                else:
                    total_spc=np.sqrt((vaxis-vsys)**2)/throw < 3
                im_mask=np.sum(mask3d, axis=0)
                # Here checking for bad fits. Places where parameter uncertainties are zero or unreasonably low
                # Probably a more elegant way to do this! 
                # For Gaussian fits: parameters are [amp, vlsr, sigma]
                for ii in np.arange( im_mask.shape[1]):
                    for jj in np.arange( im_mask.shape[0]):
                        if ((im_mask[jj,ii] == 0) or (pycube.parcube[2,jj,ii] < 3*pycube.errcube[2,jj,ii]) or
                            (pycube.errcube[1,jj,ii] == 0) or (pycube.errcube[1,jj,ii] > 0.2)):
                            # Find vlsr of nearby fits
                            vmap_loc = vmap[jj-20:jj+20,ii-20:ii+20]
                            vmean_loc = np.mean(vmap_loc[vmap_loc != 0])*u.km/u.s
                            if np.isfinite(vmean_loc):
                                vshift = np.int((vmean - vmean_loc)/dv)
                                # Shift total_spc to match local vlsr
                                total_spc_shift = np.roll(total_spc,vshift)
                            else:
                                # Try bigger window
                                vmap_loc = vmap[jj-50:jj+50,ii-50:ii+50]
                                vmean_loc = np.mean(vmap_loc[vmap_loc != 0])*u.km/u.s
                                if np.isfinite(vmean_loc):
                                    vshift = np.int((vmean - vmean_loc)/dv)
                                    # Shift total_spc to match local vlsr
                                    total_spc_shift = np.roll(total_spc,vshift)
                                else:
                                    # Use original window if no good fits nearby
                                    total_spc_shift = total_spc
                            mask3d[:,jj,ii] = total_spc_shift
                n_chan=np.sum(mask3d, axis=0)
                # create masked cube
                cube2 = cube.with_mask(mask3d)
                cube3 = cube.with_mask(~mask3d)
                #
                if save_masked:
                    cube2.write( file_temp, overwrite=True)
                # calculate moment map
                moment_0 = cube2.moment(axis=0)
                moment_0.write( file_out, overwrite=True)
                cube3.allow_huge_operations=True
                rms=cube3.std(axis=0)
                rms.write( file_rms, overwrite=True)
                mom_0_rms=rms * dv * np.sqrt(n_chan)
                mom_0_rms.write( file_rms_mom, overwrite=True)
            except IOError:
                warnings.warn("File not found: {0}".format(fit_model_file))

def run_moment_rest_all(regions=None,file_extension='all_rebase3',threshold=0.0125,save_masked=False,trim_edge=True,release='all'):
    """
    Update the moment calculations for all regions for non-NH3 lines
    Get list of regions - run from images/ directory
    Assume directories correspond to regions to be imaged
    Update - use catalog?
    """
    if file_extension is None:
        file_extension = '_'+release

    if regions is None:
        RegionCatalog = catalogs.GenerateRegions(release=release)
    else:
        RegionCatalog = catalogs.GenerateRegions(release=release)
        keep = [idx for idx, row in enumerate(RegionCatalog) if row['Region name'] in regions]
        RegionCatalog = RegionCatalog[keep]

    for ThisRegion in RegionCatalog:
        region_name=ThisRegion['Region name']
        print('Now {0}'.format(region_name))
        vsys = ThisRegion['VAVG']*u.km/u.s
        throw = 2.*u.km/u.s + ThisRegion['VRANGE']*u.km/u.s/2.
        update_rest_moment0_2(region_name=region_name,
                              file_extension=file_extension,
                              threshold=0.0125, save_masked=save_masked,
                              trim_edge=trim_edge,vsys=vsys,throw=throw)


def run_plot_fit_all():
    """
    Run the functions for fitting the NH3 line profile and plotting the 
    results for all regions
    """
    cubefit(region='L1455', blorder=1, do_plot=True, snr_min=3, multicore=1, 
            vmax=6.7, vmin=3.5, file_extension='base_DR1')
    plot_cubefit(region='L1455', distance=250*u.pc, dvmin=0.05, dvmax=0.45, 
                 vcmin=4.5, vcmax=6.0, file_extension='base_DR1')

    cubefit(region='NGC1333', blorder=1, do_plot=True, snr_min=3, multicore=1,
            vmax=9.5, vmin=4.2, file_extension='DR1_rebase3')
    plot_cubefit(region='NGC1333', distance=250*u.pc, dvmin=0.05, dvmax=0.6, 
                 vcmin=6.4, vcmax=9.3, file_extension='DR1_rebase3')

    cubefit(region='B18', vmin=4.5, vmax=7.5, do_plot=False, snr_min=3.0, 
            multicore=1, file_extension='DR1_rebase3')
    plot_cubefit(region='B18', distance=137*u.pc, dvmin=0.05, dvmax=0.3, 
                 vcmin=5.7, vcmax=6.7, file_extension='DR1_rebase3')

    cubefit(region='L1688', vmin=2.5, vmax=5.5, do_plot=False, snr_min=3.0, 
            multicore=1, file_extension='DR1_rebase3')
    plot_cubefit(region='L1688', distance=120*u.pc, dvmin=0.05, dvmax=0.7, 
                 vcmin=2.7, vcmax=4.8, file_extension='DR1_rebase3')

    cubefit(region='OrionA', vmin=5.6, vmax=13.7, do_plot=False, snr_min=3.0, 
            multicore=1, file_extension='DR1_rebase3')
    plot_cubefit(region='OrionA', distance=450*u.pc, dvmin=0.05, dvmax=0.7, 
                 vcmin=5.7, vcmax=12.7, file_extension='DR1_rebase3')


def trim_cubes(region_name='OrionA',file_extension='DR1_rebase3',blorder=1,propertyMaps=True):
    '''
    Trim cube edges for DR1 release data. Wasn't incorporated into baseline fitting
    but should be for future baselined data. 
    Since trim looks for edges in the data, CANNOT use on already flagged data.
    Use updated moment map as mask for other files
    '''
    if file_extension:
        root = file_extension
    else:
        root = '{0}'.format(blorder)

    # Cubes:
    line_list=['NH3_11','NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    #line_list = ['NH3_11','NH3_22']
    for line in line_list:
        # Moment first to create mask
        moment = fits.open('{0}/{0}_{1}_{2}_mom0_QA.fits'.format(region_name,line,file_extension))
        moment_data = moment[0].data
        moment_hdr  = moment[0].header
        moment.close()
        trim_edge_cube(moment_data)
        mask = np.isfinite(moment_data)
        # Write out new moment
        fits.writeto('{0}/{0}_{1}_{2}_mom0_QA_trim.fits'.format(region_name,line,file_extension),
                     moment_data,moment_hdr,clobber=True)
        # Next cubes
        filein = '{0}/{0}_{1}_{2}.fits'.format(region_name,line,file_extension)
        # trim_edge_cube doesn't work on a spectral cube object. 
        cube = SpectralCube.read(filein)
        cube2 = cube.with_mask(mask)
        cube2.write('{0}/{0}_{1}_{2}_trim.fits'.format(region_name,line,file_extension),overwrite=True)
        # And rms. Note that QA rms file has interior masked regions for Orion A, others?
        rms = fits.open('{0}/{0}_{1}_{2}_rms_QA.fits'.format(region_name,line,file_extension))
        rms_data = rms[0].data
        rms_hdr  = rms[0].header
        rms.close()
        trim_edge_cube(rms_data)
        fits.writeto('{0}/{0}_{1}_{2}_rms_QA_trim.fits'.format(region_name,line,file_extension),
                     rms_data,rms_hdr,clobber=True)

    # Use NH3 (1,1) moment map as mask for property map
    if propertyMaps:
        fit_file = '{0}/{0}_parameter_maps_{1}.fits'.format(region_name,file_extension)
        moment = fits.open('{0}/{0}_NH3_11_{1}_mom0_QA_trim.fits'.format(region_name,file_extension))
        moment_data = moment[0].data
        propMap = fits.open(fit_file)
        propMap_data = propMap[0].data
        propMap_hdr  = propMap[0].header
        propMap.close()    
        for plane_i in range(len(propMap_data)):
            data_i = propMap_data[plane_i,:,:]
            data_i[~np.isfinite(moment_data)] = np.nan
            propMap_data[plane_i,:,:] = data_i
        fits.writeto('{0}/{0}_parameter_maps_{1}_trim.fits'.format(region_name,file_extension),
                     propMap_data,propMap_hdr,clobber=True)
    

def _add_plot_text( fig, region, blorder, distance):
    # set nan color and beam size color
    fig.set_nan_color('0.8')
    fig.add_beam()
    fig.beam.set_color('#E31A1C')
    # Scale bar
    ang_sep = (0.1*u.pc/distance)*u.rad
    fig.add_scalebar(ang_sep)
    fig.scalebar.set(color='black')
    fig.scalebar.set_label('0.1 pc')
    # Add region and tick labels parameters
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm')
    fig.add_label(0.05,0.90, region, relative=True, color='black')
    fig.ticks.set_color('black')
    fig.ticks.set_minor_frequency(4)

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

def flag_all_data(region='OrionA',blorder='1', file_extension='DR1_rebase3'):
    """
    Flag cubefit results based on S/N in integrated intensity.
    Also flag poorly constrained fits (where Tk, Tex hit minimum values)
    Outputs individual .fits files for parameters and uncertainties: Tkin, Tex, Vc, sigmaV, NNH3
    Remove moment 0 flagging as update_NH3_moment + edge masking is better
    Update (16/9/13): Include masking of 'island pixels' based on S/N in NH3 (1,1) moment map
                      Set flagged values to zero rather than NaN
                      Output masked cube file
                      Cleaned up

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

    flagSN22 = 3.0
    flagSN11 = 3.0

    if file_extension:
        root = file_extension
    else:
        root = 'base{0}'.format(blorder)

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
    rms22hdu = fits.open("{0}/{0}_NH3_22_{1}_rms.fits".format(region,root))
    rms22data = rms22hdu[0].data
    rms22hdu.close()
    m0_22 = fits.open("{0}/{0}_NH3_22_{1}_Tpeak.fits".format(region,root))
    m0_22data = m0_22[0].data
    hd22 = m0_22[0].header
    m0_22.close()

    sn11 = m0_11data/rms11data
    sn11HDU = fits.PrimaryHDU(sn11,hd11)
    sn22 = m0_22data/rms22data

    # Get binary mask
    selem=np.array([[0,1,0],[1,1,1],[0,1,0]])
    LowestContour=3
    pixel_mask = mask_binary(sn11HDU,LowestContour=LowestContour,selem=selem)

    # Get Tex, Tk files for mask
    tex  = np.copy(cube[1,:,:])
    etex = np.copy(cube[7,:,:])
    tk   = np.copy(cube[0,:,:])
    etk  = np.copy(cube[6,:,:])

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
    parMask = ((sn22 >= flagSN22) & (tex >flagMinTex) & (param >flagMinTk) & (param <1e3))
    param = param * pixel_mask * parMask
    eparam = eparam * pixel_mask * parMask
    file_out="{0}/parameterMaps/{0}_Tkin_{1}_flag.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    # Write out uncertainties
    file_out="{0}/parameterMaps/{0}_eTkin_{1}_flag.fits".format(region,root)
    fits.writeto(file_out, eparam, hd, clobber=True)
    # Update cube
    cube[0,:,:] = param
    cube[6,:,:] = eparam

    #Tex
    hd['BUNIT']='K'
    param=cube[1,:,:]
    eparam=cube[7,:,:]
    parMask = ((sn11>=flagSN11) & (tex > flagMinTex) & (tk > flagMinTk))
    param = param * pixel_mask * parMask
    eparam = eparam * pixel_mask * parMask
    file_out="{0}/parameterMaps/{0}_Tex_{1}_flag.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    # Write out uncertainties
    file_out="{0}/parameterMaps/{0}_eTex_{1}_flag.fits".format(region,root)
    fits.writeto(file_out, eparam, hd, clobber=True)
    # Update cube
    cube[1,:,:] = param
    cube[7,:,:] = eparam

    # N_NH3
    hd['BUNIT']='cm-2'
    param=cube[2,:,:]
    eparam=cube[8,:,:]
    parMask = ((sn22 >= flagSN22) & (tex > flagMinTex) & (tk > flagMinTk))
    param = param * pixel_mask * parMask
    eparam = eparam * pixel_mask * parMask
    file_out="{0}/parameterMaps/{0}_N_NH3_{1}_flag.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    # Write out uncertainties
    file_out="{0}/parameterMaps/{0}_eN_NH3_{1}_flag.fits".format(region,root)
    fits.writeto(file_out, eparam, hd, clobber=True)    # sigma
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
    file_out="{0}/parameterMaps/{0}_Sigma_{1}_flag.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    # Write out uncertainties
    file_out="{0}/parameterMaps/{0}_eSigma_{1}_flag.fits".format(region,root)
    fits.writeto(file_out, eparam, hd, clobber=True)
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
    file_out="{0}/parameterMaps/{0}_Vlsr_{1}_flag.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    # Write out uncertainties
    file_out="{0}/parameterMaps/{0}_eVlsr_{1}_flag.fits".format(region,root)
    fits.writeto(file_out, eparam, hd, clobber=True)
    # Update cube
    cube[4,:,:] = param
    cube[10,:,:] = eparam

    # Write out new cube
    cube_out="{0}/{0}_parameter_maps_{1}_flag.fits".format(region,root)
    fits.writeto(cube_out,cube,hd_cube,clobber=True)

def plot_cubefit(region='NGC1333', blorder=1, distance=145*u.pc, dvmin=0.05, 
                 dvmax=None, vcmin=None, vcmax=None, file_extension=None):
    """
    Plot fit parameters of NH3 map for the region. It masks poorly 
    constrained fits. Extra parameters to improve images of centroid velocity
    and velocity dispersion.
    It creates 4 PDF files: Tkin, Tex, Vlsr, and sigma_v.
    
    Parameters
    ----------
    region : str
        Name of region to reduce
    blorder : int
        order of baseline removed
    distance : astropy.units
        distance to the observed region. Used for scalebar. Default 145*u.pc.
    dvmin : numpy.float
        Minimum velocity dispersion to plot, in km/s. Default 0.05
    dvmax : numpy.float
        Maximum velocity dispersion to plot, in km/s.
    vcmin : numpy.float
        Minimum centroid velocity to plot, in km/s.
    vcmax : numpy.float
        Maximum centroid velocity to plot, in km/s.
    """
    import matplotlib
    matplotlib.use('Agg')
    import aplpy
    
    if file_extension:
        root = file_extension
    else:
        root = '{0}'.format(blorder)

    data_file ="{0}_parameter_maps_{1}.fits".format(region,root)
    w11_file='{0}/{0}_NH3_11_base{1}_mom0.fits'.format(region,root)

    hdu   =fits.open(data_file)
    header=hdu[0].header
    data  =hdu[0].data
    hdu.close()
    rm_key=['NAXIS3','CRPIX3','CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3']
    for key_i in rm_key:
        header.remove(key_i)
    header['NAXIS'] = 2
    header['WCSAXES'] = 2
    # 
    # Create masked Centroid velocity map
    vc = data[4,:,:] 
    evc = data[10,:,:] 
    vc[ evc > 0.1] = np.nan
    vc[ vc == 0.0] = np.nan
    hdu_vc = fits.PrimaryHDU(vc, header)
    # Create masked velocity dispersion map
    dv = data[3,:,:] 
    edv = data[9,:,:] 
    dv[ edv > 0.2*dv] = np.nan
    #dv[ edv > 0.1] = np.nan
    dv[ evc > 0.1] = np.nan
    dv[ dv == 0.0] = np.nan
    hdu_dv = fits.PrimaryHDU(dv, header)
    # Create masked Kinetic Temperature map
    tk = data[0,:,:] 
    etk = data[6,:,:] 
    tk[ etk > 3] = np.nan
    tk[ tk == 0.0] = np.nan
    hdu_tk = fits.PrimaryHDU(tk, header)
    # Create masked Excitation Temperature map
    tex = data[1,:,:] 
    etex = data[7,:,:] 
    tex[ etex > 1.0] = np.nan
    tex[ tex == 0.0] = np.nan
    hdu_tex = fits.PrimaryHDU(tex, header)
    # Create masked N(NH3) map
    nnh3 = data[2,:,:]
    ennh3 = data[8,:,:]
    nnh3[ ennh3 > 0.3*nnh3] = np.nan
    nnh3[ nnh3 == 0] = np.nan
    hdu_nnh3 = fits.PrimaryHDU(nnh3, header)
    
    c_levs=np.arange(0.3,5,0.5)
    #
    # Centroid velocity
    #
    color_table='RdYlBu_r'
    fig0=aplpy.FITSFigure(hdu_vc, hdu=0)
    fig0.show_colorscale( cmap=color_table, vmin=vcmin, vmax=vcmax)
    fig0.show_contour(w11_file, colors='black', linewidths=0.5, levels=c_levs,
                      zorder=34)
    _add_plot_text( fig0, region, blorder, distance)
    # Colorbar 
    fig0.add_colorbar()
    fig0.colorbar.set_location('right')
    fig0.colorbar.set_axis_label_text('V$_{LSR}$ (km s$^{-1}$)')
    # Save file
    fig0.save("{0}_Vc.pdf".format(region),dpi=300)
    fig0.close()
    #
    # Sigma
    # 
    color_table='Blues'
    fig0=aplpy.FITSFigure(hdu_dv, hdu=0)
    fig0.show_colorscale( cmap=color_table, vmin=dvmin, vmax=dvmax)
    fig0.show_contour(w11_file, colors='black', linewidths=0.5, levels=c_levs,
                      zorder=34)
    _add_plot_text( fig0, region, blorder, distance)
    # Colorbar 
    fig0.add_colorbar()
    fig0.colorbar.set_location('right')
    fig0.colorbar.set_axis_label_text('$\sigma_{v}$ (km s$^{-1}$)')
    # Save file
    fig0.save("{0}_sigmaV.pdf".format(region),dpi=300)
    fig0.close()
    #
    # Tkin
    # 
    color_table='YlOrBr'
    fig0=aplpy.FITSFigure(hdu_tk, hdu=0)
    fig0.show_colorscale( cmap=color_table)
    fig0.show_contour(w11_file, colors='black', linewidths=0.5, levels=c_levs,
                      zorder=34)
    _add_plot_text( fig0, region, blorder, distance)
    # Colorbar 
    fig0.add_colorbar()
    fig0.colorbar.set_location('right')
    fig0.colorbar.set_axis_label_text('T$_{kin}$ (K)')
    # Save file
    fig0.save("{0}_Tkin.pdf".format(region),dpi=300)
    fig0.close()
    #
    # Tex
    # 
    color_table='YlOrBr'
    fig0=aplpy.FITSFigure(hdu_tex, hdu=0)
    fig0.show_colorscale( cmap=color_table)
    fig0.show_contour(w11_file, colors='black', linewidths=0.5, levels=c_levs,
                      zorder=34)
    _add_plot_text( fig0, region, blorder, distance)
    # Colorbar 
    fig0.add_colorbar()
    fig0.colorbar.set_location('right')
    fig0.colorbar.set_axis_label_text('T$_{ex}$ (K)')
    # Save file
    fig0.save("{0}_Tex.pdf".format(region),dpi=300)
    fig0.close()
    #
    # N(NH3)
    # 
    color_table='viridis'
    fig0=aplpy.FITSFigure(hdu_nnh3, hdu=0)
    fig0.show_colorscale( cmap=color_table)
    fig0.show_contour(w11_file, colors='black', linewidths=0.5, levels=c_levs,
                      zorder=34)
    _add_plot_text( fig0, region, blorder, distance)
    # Colorbar 
    fig0.add_colorbar()
    fig0.colorbar.set_location('right')
    fig0.colorbar.set_axis_label_text('N(NH$_3$) (cm$^{-2}$)')
    # Save file
    fig0.save("{0}_NNH3.pdf".format(region),dpi=300)
    fig0.close()

def plot_all_flagged(region='OrionA', blorder=1, distance=450.*u.pc, 
                     version='v1',dvmin=0, file_extension='DR1_rebase3',
                     dvmax=None, vcmin=None, vcmax=None):
    """
    Plot from flagged fits files rather than from cubefit output multi-HDU 
    file. Side-by-side plots for v_lsr, sigma
    Parameters
    ----------
    region : str
        Name of region to reduce
    blorder : int
        order of baseline removed
    distance : astropy.units
        distance to the observed region. Used for scalebar. Default 145*u.pc.
    version : str
        flagged data will have 'v1', 'v2', etc. 
    dvmin : numpy.float
        Minimum velocity dispersion to plot, in km/s. No default. If flagged
        well, shouldn't need this or other parameters below.
    dvmax : numpy.float
        Maximum velocity dispersion to plot, in km/s. No default.
    vcmin : numpy.float
        Minimum centroid velocity to plot, in km/s. No default.
    vcmax : numpy.float
        Maximum centroid velocity to plot, in km/s. No default.
    """
    if file_extension:
        root = file_extension
    else:
        root = '{0}'.format(blorder)

    w11_file = "{0}/{0}_NH3_11_{1}_mom0_flag.fits".format(region,root)
    c_levs=np.arange(0.3,5,0.5)
    #
    # Centroid velocity
    #
    # Masked file
    color_table='RdYlBu_r'
    dataFile = "{0}/parameterMaps/{0}_Vlsr_{1}_flag.fits".format(region,root)
    fig1=aplpy.FITSFigure(dataFile)
    fig1.show_colorscale( cmap=color_table, vmin=vcmin, vmax=vcmax)
    fig1.show_contour(w11_file, colors='black', linewidths=0.5, levels=c_levs,
                      zorder=34)
    fig1.tick_labels.set_style('colons')
    fig1.tick_labels.set_yformat('dd:mm')    
    fig1.tick_labels.set_xformat('hh:mm')
    # Colorbar
    fig1.add_colorbar()
    fig1.colorbar.set_location('right')
    fig1.colorbar.set_axis_label_text('V$_\mathrm{LSR}$ (km s$^{-1}$)')    
    # Save file
    fig1.save("{0}/parameterMaps/{0}_Vlsr_{1}_flag.pdf".format(region,root),dpi=300)
    fig1.close()
    #
    # Sigma
    #
    # Masked file
    dataFile = "{0}/parameterMaps/{0}_Sigma_{1}_flag.fits".format(region,root)
    fig1=aplpy.FITSFigure(dataFile)
    fig1.show_colorscale( cmap=color_table, vmin=dvmin, vmax=dvmax)
    fig1.show_contour(w11_file, colors='black', linewidths=0.5, levels=c_levs,
                      zorder=34)
    fig1.tick_labels.set_style('colons')
    fig1.tick_labels.set_yformat('dd:mm')
    fig1.tick_labels.set_xformat('hh:mm')    
    # Colorbar 
    fig1.add_colorbar()
    fig1.colorbar.set_location('right')
    fig1.colorbar.set_axis_label_text('$\sigma_\mathrm{v}$ (km s$^{-1}$)')
    # Save file
    fig1.save("{0}/{0}_Sigma_{1}_flag.pdf".format(region,root),dpi=300)
    fig1.close()
    #
    # Tkin
    # 
    dataFile = "{0}/parameterMaps/{0}_Tkin_{1}_flag.fits".format(region,root)
    color_table='YlOrBr'
    fig0=aplpy.FITSFigure(dataFile)
    fig0.show_colorscale( cmap=color_table)
    fig0.show_contour(w11_file, colors='black', linewidths=0.5, levels=c_levs,
                      zorder=34)
    _add_plot_text( fig0, region, blorder, distance)
    # Colorbar 
    fig0.add_colorbar()
    fig0.colorbar.set_location('right')
    fig0.colorbar.set_axis_label_text('T$_\mathrm{kin}$ (K)')
    # Save file
    fig0.save("{0}/parameterMaps/{0}_Tkin_{1}_flag.pdf".format(region,root),dpi=300)
    fig0.close()
    #
    # Tex
    # 
    dataFile = "{0}/parameterMaps/{0}_Tex_{1}_flag.fits".format(region,root)
    color_table='YlOrBr'
    fig0=aplpy.FITSFigure(dataFile)
    fig0.show_colorscale( cmap=color_table)
    fig0.show_contour(w11_file, colors='black', linewidths=0.5, levels=c_levs,
                      zorder=34)
    _add_plot_text( fig0, region, blorder, distance)
    # Colorbar 
    fig0.add_colorbar()
    fig0.colorbar.set_location('right')
    fig0.colorbar.set_axis_label_text('T$_\mathrm{ex}$ (K)')
    # Save file
    fig0.save("{0}/parameterMaps/{0}_Tex_{1}_flag.pdf".format(region,root),dpi=300)
    fig0.close()
    #
    # N(NH3)
    # 
    dataFile = "{0}/parameterMaps/{0}_N_NH3_{1}_flag.fits".format(region,root)
    color_table='YlOrBr'
    fig0=aplpy.FITSFigure(dataFile)
    fig0.show_colorscale( cmap=color_table)
    fig0.show_contour(w11_file, colors='black', linewidths=0.5, levels=c_levs,
                      zorder=34)
    _add_plot_text( fig0, region, blorder, distance)
    # Colorbar 
    fig0.add_colorbar()
    fig0.colorbar.set_location('right')
    fig0.colorbar.set_axis_label_text('N(NH$_3$) (cm$^{-2}$)')
    # Save file
    fig0.save("{0}/parameterMaps/{0}_NNH3_{1}_flag.pdf".format(region,root),dpi=300)
    fig0.close()

def update_cubefit(region='NGC1333', blorder=1, file_extension=None):
    """
    Updates the fit parameters storage format from cube (v0, one channel per 
    parameter) into a set of files (v1, one FITS per parameter). 
    """
    if file_extension:
        root = file_extension
    else:
        root = 'base{0}'.format(blorder)

    hdu=fits.open("{0}/{0}_parameter_maps_{1}.fits".format(region,root))
    hd=hdu[0].header
    cube=hdu[0].data
    hdu.close()

    rm_key=['NAXIS3','CRPIX3','CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3']
    for key_i in rm_key:
        hd.remove(key_i)
    hd['NAXIS']= 2
    hd['WCSAXES']= 2
    # Tkin
    hd['BUNIT']='K'
    param=cube[0,:,:]
    file_out="{0}/{0}_Tkin_{1}.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    #Tex
    hd['BUNIT']='K'
    param=cube[1,:,:]
    file_out="{0}/{0}_Tex_{1}.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    # N_NH3
    hd['BUNIT']='cm-2'
    param=cube[2,:,:]
    file_out="{0}/{0}_N_NH3_{1}.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    # sigma
    hd['BUNIT']='km/s'
    param=cube[3,:,:]
    file_out="{0}/{0}_Sigma_{1}.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    # Vlsr
    hd['BUNIT']='km/s'
    param=cube[4,:,:]
    file_out="{0}/{0}_Vlsr_{1}.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    # Fortho
    hd['BUNIT']=''
    param=cube[5,:,:]
    file_out="{0}/{0}_Fortho_{1}.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    # eTkin
    hd['BUNIT']='K'
    param=cube[6,:,:]
    file_out="{0}/{0}_eTkin_{1}.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    # eTex
    hd['BUNIT']='K'
    param=cube[7,:,:]
    file_out="{0}/{0}_eTex_{1}.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    # eN_NH3
    hd['BUNIT']='cm-2'
    param=cube[8,:,:]
    file_out="{0}/{0}_eN_NH3_{1}.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    # eSigma
    hd['BUNIT']='km/s'
    param=cube[9,:,:]
    file_out="{0}/{0}_eSigma_{1}.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    # eVlsr
    hd['BUNIT']='km/s'
    param=cube[10,:,:]
    file_out="{0}/{0}_eVlsr_{1}.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)
    # eFortho
    hd['BUNIT']=''
    param=cube[11,:,:]
    file_out="{0}/{0}_eFortho_{1}.fits".format(region,root)
    fits.writeto(file_out, param, hd, clobber=True)

def default_masking(snr,snr_min=5.0):
    planemask = (snr>snr_min) 
    planemask = remove_small_objects(planemask,min_size=40)
    planemask = opening(planemask,disk(1))
    return(planemask)


def cubefit(region='NGC1333', blorder=1, vmin=5, vmax=15, do_plot=False, 
            snr_min=5.0, multicore=1, file_extension=None, mask_function = None, gauss_fit=False):
    """
    Fit NH3(1,1) and (2,2) cubes for the requested region. 
    It fits all pixels with SNR larger than requested. 
    Initial guess is based on moment maps and neighboring pixels. 
    The fitting can be done in parallel mode using several cores, 
    however, this is dangerous for large regions, where using a 
    good initial guess is important. 
    It stores the result in a FITS cube. 

    TODO:
    -Improve initial guess
    
    Parameters
    ----------
    region : str
        Name of region to reduce
    blorder : int
        order of baseline removed
    vmin : numpy.float
        Minimum centroid velocity to plot, in km/s.
    vmax : numpy.float
        Maximum centroid velocity to plot, in km/s.
    do_plot : bool
        If True, then a map of the region to map is shown.
    snr_min : numpy.float
        Minimum signal to noise ratio of the spectrum to be fitted.
    multicore : int
        Numbers of cores to use for parallel processing.
    file_extension : str
        File extension of the input maps. Default is 'base#' where # is the 
        blorder parameter above.
    mask_function : fun
        function to create a custom made mask for analysis. Defaults to using 
        `default_masking`
    """
    if file_extension:
        root = file_extension
    else:
        # root = 'base{0}'.format(blorder)
        root = 'all'

    OneOneIntegrated = '{0}/{0}_NH3_11_{1}_mom0.fits'.format(region,root)
    OneOneFile = '{0}/{0}_NH3_11_{1}.fits'.format(region,root)
    RMSFile = '{0}/{0}_NH3_11_{1}_rms.fits'.format(region,root)
    TwoTwoFile = '{0}/{0}_NH3_22_{1}.fits'.format(region,root)
    ThreeThreeFile = '{0}/{0}_NH3_33_{1}.fits'.format(region,root)
        
    cube11sc = SpectralCube.read(OneOneFile)
    cube22sc = SpectralCube.read(TwoTwoFile)
    errmap11 = fits.getdata(RMSFile)
    rms = np.nanmedian(errmap11)

    snr = cube11sc.filled_data[:].value/errmap11
    peaksnr = np.max(snr,axis=0)
    if mask_function is None:
        planemask = default_masking(peaksnr,snr_min = snr_min)
    else:
        planemask = mask_function(peaksnr,snr_min = snr_min)
    
    #planemask = (peaksnr>20) * (errmap11 < 0.2)

    mask = (snr>3)*planemask
    maskcube = cube11sc.with_mask(mask.astype(bool))
    maskcube = maskcube.with_spectral_unit(u.km/u.s,velocity_convention='radio')
    slab = maskcube.spectral_slab( vmax*u.km/u.s, vmin*u.km/u.s)
    w11=slab.moment( order=0, axis=0).value
    peakloc = np.nanargmax(w11)
    ymax,xmax = np.unravel_index(peakloc,w11.shape)
    moment1 = slab.moment( order=1, axis=0).value
    moment2 = (slab.moment( order=2, axis=0).value)**0.5
    moment2[np.isnan(moment2)]=0.2
    moment2[moment2<0.2]=0.2
    cube11 = pyspeckit.Cube(OneOneFile,maskmap=planemask)
    cube11.unit="K"
    cube22 = pyspeckit.Cube(TwoTwoFile,maskmap=planemask)
    cube22.unit="K"
    #cube33 = pyspeckit.Cube(ThreeThreeFile,maskmap=planemask)
    #cube33.unit="K" # removed as long as we're not modeling OPR
    cubes = pyspeckit.CubeStack([cube11,cube22],maskmap=planemask)
    cubes.unit="K"
    guesses = np.zeros((6,)+cubes.cube.shape[1:])
    moment1[moment1<vmin] = vmin+0.2
    moment1[moment1>vmax] = vmax-0.2
    guesses[0,:,:] = 12                    # Kinetic temperature 
    guesses[1,:,:] = 3                     # Excitation  Temp
    guesses[2,:,:] = 14.5                  # log(column)
    guesses[3,:,:] = moment2  # Line width / 5 (the NH3 moment overestimates linewidth)               
    guesses[4,:,:] = moment1  # Line centroid              
    guesses[5,:,:] = 0.0                   # F(ortho) - ortho NH3 fraction (fixed)
    if do_plot:
        import matplotlib.pyplot as plt
        plt.imshow( w11, origin='lower',interpolation='nearest')
        plt.show()
    F=False
    T=True
    
    if not 'cold_ammonia' in cubes.specfit.Registry.multifitters:
        cubes.specfit.Registry.add_fitter('cold_ammonia',ammonia.cold_ammonia_model(),6)
        
    print('start fit')
    cubes.fiteach(fittype='cold_ammonia',  guesses=guesses,
                  integral=False, verbose_level=3, 
                  fixed=[F,F,F,F,F,T], signal_cut=2,
                  limitedmax=[F,F,T,F,T,T],
                  maxpars=[0,0,17.0,0,vmax,1],
                  limitedmin=[T,T,T,T,T,T],
                  minpars=[5,2.8,12.0,0.04,vmin,0],
                  start_from_point=(xmax,ymax),
                  use_neighbor_as_guess=False, 
                  position_order = 1/peaksnr,
                  errmap=errmap11, multicore=multicore)

    fitcubefile = fits.PrimaryHDU(data=np.concatenate([cubes.parcube,cubes.errcube]), header=cubes.header)
    fitcubefile.header.set('PLANE1','TKIN')
    fitcubefile.header.set('PLANE2','TEX')
    fitcubefile.header.set('PLANE3','COLUMN')
    fitcubefile.header.set('PLANE4','SIGMA')
    fitcubefile.header.set('PLANE5','VELOCITY')
    fitcubefile.header.set('PLANE6','FORTHO')
    fitcubefile.header.set('PLANE7','eTKIN')
    fitcubefile.header.set('PLANE8','eTEX')
    fitcubefile.header.set('PLANE9','eCOLUMN')
    fitcubefile.header.set('PLANE10','eSIGMA')
    fitcubefile.header.set('PLANE11','eVELOCITY')
    fitcubefile.header.set('PLANE12','eFORTHO')
    fitcubefile.header.set('CDELT3',1)
    fitcubefile.header.set('CTYPE3','FITPAR')
    fitcubefile.header.set('CRVAL3',0)
    fitcubefile.header.set('CRPIX3',1)
    fitcubefile.writeto("{0}/{0}_parameter_maps_{1}.fits".format(region,root),clobber=True)

    if gauss_fit==True:
	molecules = ['C2S', 'HC7N_22_21', 'HC7N_21_20', 'HC5N']
	for i in molecules:
        	gauss_fitter(region=region, mol=i, vmin=vmin, vmax=vmax, snr_min=snr_min, multicore=multicore, file_extension=file_extension)

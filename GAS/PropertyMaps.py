import pyspeckit
import astropy.io.fits as fits
import numpy as np
import os
from spectral_cube import SpectralCube
import signal_id
from radio_beam import Beam
import astropy.constants as con
import astropy.units as u
from skimage.morphology import remove_small_objects,closing,disk,opening

def run_plot_fit_all():
    """
    Run the functions for fitting the NH3 line profile and plotting the 
    results for all regions
    """
    cubefit(region='L1455', blorder=1, do_plot=True, snr_min=3, multicore=1, 
            vmax=6.7, vmin=3.5)
    plot_cubefit(region='L1455', distance=250*u.pc, dvmin=0.05, dvmax=0.45, 
                 vcmin=4.5, vcmax=6.0)

    cubefit(region='NGC1333', blorder=1, do_plot=True, snr_min=3, multicore=1,
            vmax=9.5, vmin=4.2)
    plot_cubefit(region='NGC1333', distance=250*u.pc, dvmin=0.05, dvmax=0.6, 
                 vcmin=6.4, vcmax=9.3)

    cubefit(region='B18', vmin=4.5, vmax=7.5, do_plot=False, snr_min=3.0, 
            multicore=1)
    plot_cubefit(region='B18', distance=137*u.pc, dvmin=0.05, dvmax=0.3, 
                 vcmin=5.7, vcmax=6.7)

    cubefit(region='L1688', vmin=5.5, vmax=10.5, do_plot=False, snr_min=3.0, 
            multicore=1)
    plot_cubefit(region='L1688', distance=120*u.pc, dvmin=0.05, dvmax=0.7, 
                 vcmin=2.7, vcmax=4.8)

    cubefit(region='OrionA', vmin=5.6, vmax=13.7, do_plot=False, snr_min=3.0, 
            multicore=1)
    plot_cubefit(region='OrionA', distance=450*u.pc, dvmin=0.05, dvmax=0.7, 
                 vcmin=5.7, vcmax=12.7)



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

def plot_cubefit(region='NGC1333', blorder=1, distance=145*u.pc, dvmin=0.05, 
                 dvmax=None, vcmin=None, vcmax=None):
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

    data_file = "{0}_parameter_maps.fits".format(region)
    w11_file='{0}/{0}_NH3_11_base{1}_mom0.fits'.format(region,blorder)

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
    dv[ edv > 0.1] = np.nan
    dv[ dv == 0.0] = np.nan
    hdu_dv = fits.PrimaryHDU(dv, header)
    # Create masked Kinetic Temperature map
    tk = data[0,:,:] 
    etk = data[6,:,:] 
    tk[ etk > 0.5] = np.nan
    tk[ tk == 0.0] = np.nan
    hdu_tk = fits.PrimaryHDU(tk, header)
    # Create masked Excitation Temperature map
    tex = data[1,:,:] 
    etex = data[7,:,:] 
    tex[ etex > 1.0] = np.nan
    tex[ tex == 0.0] = np.nan
    hdu_tex = fits.PrimaryHDU(tex, header)
    
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
    # Tkin
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



def update_cubefit(region='NGC1333'):
    """
    Updates the fit parameters storage format from cube (v0, one channel per 
    parameter) into a set of files (v1, one FITS per parameter). 
    """
    hdu=fits.open("{0}_parameter_maps.fits".format(region))
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
    file_out="{0}_Tkin_v1.fits".format(region)
    fits.writeto(file_out, param, hd, clobber=True)
    #Tex
    hd['BUNIT']='K'
    param=cube[1,:,:]
    file_out="{0}_Tex_v1.fits".format(region)
    fits.writeto(file_out, param, hd, clobber=True)
    # N_NH3
    hd['BUNIT']='cm-2'
    param=cube[2,:,:]
    file_out="{0}_N_NH3_v1.fits".format(region)
    fits.writeto(file_out, param, hd, clobber=True)
    # sigma
    hd['BUNIT']='km/s'
    param=cube[3,:,:]
    file_out="{0}_Sigma_v1.fits".format(region)
    fits.writeto(file_out, param, hd, clobber=True)
    # Vlsr
    hd['BUNIT']='km/s'
    param=cube[4,:,:]
    file_out="{0}_Vlsr_v1.fits".format(region)
    fits.writeto(file_out, param, hd, clobber=True)
    # Fortho
    hd['BUNIT']=''
    param=cube[5,:,:]
    file_out="{0}_Fortho_v1.fits".format(region)
    fits.writeto(file_out, param, hd, clobber=True)
    # eTkin
    hd['BUNIT']='K'
    param=cube[6,:,:]
    file_out="{0}_eTkin_v1.fits".format(region)
    fits.writeto(file_out, param, hd, clobber=True)
    # eTex
    hd['BUNIT']='K'
    param=cube[7,:,:]
    file_out="{0}_eTex_v1.fits".format(region)
    fits.writeto(file_out, param, hd, clobber=True)
    # eN_NH3
    hd['BUNIT']='cm-2'
    param=cube[8,:,:]
    file_out="{0}_eN_NH3_v1.fits".format(region)
    fits.writeto(file_out, param, hd, clobber=True)
    # eSigma
    hd['BUNIT']='km/s'
    param=cube[9,:,:]
    file_out="{0}_eSigma_v1.fits".format(region)
    fits.writeto(file_out, param, hd, clobber=True)
    # eVlsr
    hd['BUNIT']='km/s'
    param=cube[10,:,:]
    file_out="{0}_eVlsr_v1.fits".format(region)
    fits.writeto(file_out, param, hd, clobber=True)
    # eFortho
    hd['BUNIT']=''
    param=cube[11,:,:]
    file_out="{0}_eFortho_v1.fits".format(region)
    fits.writeto(file_out, param, hd, clobber=True)

def cubefit(region='NGC1333', blorder=1, vmin=5, vmax=15, do_plot=False, 
            snr_min=5.0, multicore=1):
    """
    Fit NH3(1,1), (2,2) and (3,3) cubes for the requested region. 
    It fits all pixels with SNR larger than requested. 
    Initial guess is based on moment maps and neighboring pixels. 
    The fitting can be done in parallel mode using several cores, 
    however, this is dangerous for large regions, where using a 
    good initial guess is important. 
    It stores the result in a FITS cube. 

    TODO:
    -Store results in hdu list
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
    """

    OneOneIntegrated = '{0}/{0}_NH3_11_mom0.fits'.format(region,blorder)
    OneOneFile = '{0}/{0}_NH3_11_base{1}.fits'.format(region,blorder)
    RMSFile = '{0}/{0}_NH3_11_base{1}_rms.fits'.format(region,blorder)
    TwoTwoFile = '{0}/{0}_NH3_22_base{1}.fits'.format(region,blorder)
    ThreeThreeFile = '{0}/{0}_NH3_33_base{1}.fits'.format(region,blorder)
        
    beam11 = Beam.from_fits_header(fits.getheader(OneOneFile))
    cube11sc = SpectralCube.read(OneOneFile)
    cube22sc = SpectralCube.read(TwoTwoFile)
    errmap11 = fits.getdata(RMSFile)
    snr = cube11sc.filled_data[:].value/errmap11
    peaksnr = np.max(snr,axis=0)
    rms = np.nanmedian(errmap11)
    planemask = (peaksnr>snr_min) # *(errmap11 < 0.15)
    planemask = remove_small_objects(planemask,min_size=40)
    planemask = opening(planemask,disk(1))
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
    maskmap = w11>0.5
    cube11 = pyspeckit.Cube(OneOneFile,maskmap=planemask)
    cube11.unit="K"
    cube22 = pyspeckit.Cube(TwoTwoFile,maskmap=planemask)
    cube22.unit="K"
    cube33 = pyspeckit.Cube(ThreeThreeFile,maskmap=planemask)
    cube33.unit="K"
    cubes = pyspeckit.CubeStack([cube11,cube22,cube33],maskmap=planemask)
    cubes.unit="K"
    guesses = np.zeros((6,)+cubes.cube.shape[1:])
    moment1[moment1<vmin] = vmin+0.2
    moment1[moment1>vmax] = vmax-0.2
    guesses[0,:,:] = 12                    # Kinetic temperature 
    guesses[1,:,:] = 3                     # Excitation  Temp
    guesses[2,:,:] = 14.5                  # log(column)
    guesses[3,:,:] = moment2  # Line width / 5 (the NH3 moment overestimates linewidth)               
    guesses[4,:,:] = moment1  # Line centroid              
    guesses[5,:,:] = 0.5                   # F(ortho) - ortho NH3 fraction (fixed)
    if do_plot:
        import matplotlib.pyplot as plt
        plt.imshow( w11, origin='lower')
        plt.show()
    F=False
    T=True
    print('start fit')
    cubes.fiteach(fittype='ammonia',  guesses=guesses,
                  integral=False, verbose_level=3, 
                  fixed=[F,F,F,F,F,T], signal_cut=2,
                  limitedmax=[F,F,F,F,T,T],
                  maxpars=[0,0,0,0,vmax,1],
                  limitedmin=[T,T,F,T,T,T],
                  minpars=[5,2.8,0,0,vmin,0],
                  start_from_point=(xmax,ymax),
                  use_neighbor_as_guess=True, 
                  position_order = 1/peaksnr,
                  errmap=errmap11, multicore=multicore)

    fitcubefile = fits.PrimaryHDU(data=np.concatenate([cubes.parcube,cubes.errcube]), header=cubes.header)
    fitcubefile.header.update('PLANE1','TKIN')
    fitcubefile.header.update('PLANE2','TEX')
    fitcubefile.header.update('PLANE3','COLUMN')
    fitcubefile.header.update('PLANE4','SIGMA')
    fitcubefile.header.update('PLANE5','VELOCITY')
    fitcubefile.header.update('PLANE6','FORTHO')
    fitcubefile.header.update('PLANE7','eTKIN')
    fitcubefile.header.update('PLANE8','eTEX')
    fitcubefile.header.update('PLANE9','eCOLUMN')
    fitcubefile.header.update('PLANE10','eSIGMA')
    fitcubefile.header.update('PLANE11','eVELOCITY')
    fitcubefile.header.update('PLANE12','eFORTHO')
    fitcubefile.header.update('CDELT3',1)
    fitcubefile.header.update('CTYPE3','FITPAR')
    fitcubefile.header.update('CRVAL3',0)
    fitcubefile.header.update('CRPIX3',1)
    fitcubefile.writeto("{0}_parameter_maps.fits".format(region),clobber=True)

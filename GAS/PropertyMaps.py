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

def add_plot_text( fig, region, blorder, distance):
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

def plot_cubefit( region='NGC1333', blorder=1, distance=145*u.pc):
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


    #
    # Centroid velocity
    #
    color_table='RdYlBu_r'
    fig0=aplpy.FITSFigure(hdu_vc, hdu=0)
    fig0.show_colorscale( cmap=color_table)
    fig0.show_contour(w11_file, colors='black', linewidths=0.5, levels=np.arange(0.3,5,0.5), zorder=34)
    add_plot_text( fig0, region, blorder, distance)
    # Colorbar 
    fig0.add_colorbar()
    fig0.colorbar.set_location('right')
    fig0.colorbar.set_axis_label_text('V$_{LSR}$ (km s$^{-1}$)')
    # Save file
    fig0.save("{0}_Vc.pdf".format(region))
    #
    # Sigma
    # 
    color_table='Blues'
    fig0=aplpy.FITSFigure(hdu_dv, hdu=0)
    fig0.show_colorscale( cmap=color_table, vmin=0.05, vmax=0.4)
    fig0.show_contour(w11_file, colors='black', linewidths=0.5, levels=np.arange(0.3,5,0.5), zorder=34)
    add_plot_text( fig0, region, blorder, distance)
    # Colorbar 
    fig0.add_colorbar()
    fig0.colorbar.set_location('right')
    fig0.colorbar.set_axis_label_text('$\sigma_{v}$ (km s$^{-1}$)')
    # Save file
    fig0.save("{0}_sigmaV.pdf".format(region))
    #
    # Tkin
    # 
    color_table='YlOrBr'
    fig0=aplpy.FITSFigure(hdu_tk, hdu=0)
    fig0.show_colorscale( cmap=color_table)
    fig0.show_contour(w11_file, colors='black', linewidths=0.5, levels=np.arange(0.3,5,0.5), zorder=34)
    add_plot_text( fig0, region, blorder, distance)
    # Colorbar 
    fig0.add_colorbar()
    fig0.colorbar.set_location('right')
    fig0.colorbar.set_axis_label_text('T$_{kin}$ (K)')
    # Save file
    fig0.save("{0}_Tkin.pdf".format(region))
    #
    # Tkin
    # 
    color_table='YlOrBr'
    fig0=aplpy.FITSFigure(hdu_tex, hdu=0)
    fig0.show_colorscale( cmap=color_table)
    fig0.show_contour(w11_file, colors='black', linewidths=0.5, levels=np.arange(0.3,5,0.5), zorder=34)
    add_plot_text( fig0, region, blorder, distance)
    # Colorbar 
    fig0.add_colorbar()
    fig0.colorbar.set_location('right')
    fig0.colorbar.set_axis_label_text('T$_{ex}$ (K)')
    # Save file
    fig0.save("{0}_Tex.pdf".format(region))


def cubefit(region = 'NGC1333',blorder=1,vmin=5,vmax=15, do_plot=False, snr_min=5.0, multicore=1):

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

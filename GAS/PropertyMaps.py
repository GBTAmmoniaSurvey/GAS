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

def cubefit(region = 'NGC1333',blorder=1,vmin=5,vmax=15):

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
    planemask = (peaksnr>3.5)*(errmap11 < 0.15)
    planemask = remove_small_objects(planemask,min_size=40)
    planemask = opening(planemask,disk(1))
    #planemask = (peaksnr>20) * (errmap11 < 0.2)

    mask = (snr>3)*planemask
    maskcube = cube11sc.with_mask(mask.astype(bool))
    maskcube = maskcube.with_spectral_unit(u.km/u.s,velocity_convention='radio')
    slab = maskcube.spectral_slab(15*u.km/u.s,5*u.km/u.s)
    w11=slab.moment(0,0).value
    peakloc = np.nanargmax(w11)
    ymax,xmax = np.unravel_index(peakloc,w11.shape)
    moment1 = slab.moment(1,0).value
    moment2 = (slab.moment(2,0).value)**0.5
    moment2[np.isnan(moment2)]=0.2
    moment2[moment2<0.2]=0.2
    maskmap = w11>0.5
    cube11 = pyspeckit.Cube(OneOneFile,maskmap=planemask)
    cube11.units="K"
    cube22 = pyspeckit.Cube(TwoTwoFile,maskmap=planemask)
    cube22.units="K"
    cube33 = pyspeckit.Cube(ThreeThreeFile,maskmap=planemask)
    cube33.units="K"
    cubes = pyspeckit.CubeStack([cube11,cube22,cube33],maskmap=planemask)
    cubes.units="K"
    guesses = np.zeros((6,)+cubes.cube.shape[1:])
    moment1[moment1<vmin] = vmin+0.2
    moment1[moment1>vmax] = vmax-0.2
    guesses[0,:,:] = 12                    # Kinetic temperature 
    guesses[1,:,:] = 3                     # Excitation  Temp
    guesses[2,:,:] = 14.5                  # log(column)
    guesses[3,:,:] = moment2  # Line width / 5 (the NH3 moment overestimates linewidth)               
    guesses[4,:,:] = moment1  # Line centroid              
    guesses[5,:,:] = 0.5                   # F(ortho) - ortho NH3 fraction (fixed)
    F=False
    T=True
    print('start fit')
    cubes.fiteach(fittype='ammonia',  guesses=guesses,
                  integral=False, verbose_level=3, 
                  fixed=[F,F,F,F,F,T], signal_cut=2,
                  limitedmax=[F,F,F,F,T,T],
                  maxpars=[0.,0.,0.,0.,vmax,1.],
                  limitedmin=[T,T,F,T,T,T],
                  minpars=[5.,2.8,0,0,vmin,0],
                  start_from_point=(xmax,ymax),
                  use_neighbor_as_guess=True, 
                  position_order = 1/peaksnr,
                  errmap=errmap11)

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

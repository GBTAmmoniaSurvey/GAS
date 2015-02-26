import pyspeckit
import astropy.io.fits as fits
import numpy as np
import os
from spectral_cube import SpectralCube
import signal_id
from radio_beam import Beam
import astropy.constants as con
import astropy.units as u
OneOneIntegrated = 'NGC1333_W11.fits'
OneOneFile = 'NGC1333_ABCDEFGH_NH3_11_all_blsub.rg.fits'
TwoTwoFile = 'NGC1333_ABCDEFGH_NH3_22_all_blsub.fits'
ThreeThreeFile = 'NGC1333_ABCDEFGH_NH3_33_all_blsub.fits'


beam11 = Beam.from_fits_header(fits.getheader(OneOneFile))
cube11sc = SpectralCube.read(OneOneFile)
cube22sc = SpectralCube.read(TwoTwoFile)

noisy = signal_id.noise.Noise(cube11sc)
noisy.estimate_noise(spectral_flat=True,niter=2)
errmap11 = (noisy.spatial_norm*noisy.scale)
rms = np.nanmedian(errmap11)
mask = (noisy.snr>5)*(noisy.spatial_norm<2.5)
maskcube = cube11sc.with_mask(mask)
maskcube = maskcube.with_spectral_unit(u.km/u.s,velocity_convention='radio')
slab = maskcube#.spectral_slab(5*u.km/u.s,15*u.km/u.s)
w11=slab.moment(0,0).value
peakloc = np.nanargmax(w11)
xmax,ymax = np.unravel_index(peakloc,w11.shape)
moment1 = slab.moment(1,0).value
moment2 = (slab.moment(2,0).value)**0.5
moment2[np.isnan(moment2)]=0.2
moment2[moment2<0.2]=0.2

maskmap = w11>0.5
cube11 = pyspeckit.Cube(OneOneFile,maskmap=maskmap)
cube11.units="K"
cube22 = pyspeckit.Cube(TwoTwoFile,maskmap=maskmap)
cube22.units="K"
cube33 = pyspeckit.Cube(ThreeThreeFile,maskmap=maskmap)
cube33.units="K"
cubes = pyspeckit.CubeStack([cube11,cube22,cube33],maskmap=maskmap)
cubes.units="K"

guesses = np.zeros((6,)+cubes.cube.shape[1:])
guesses[0,:,:] = 12                    # Kinetic temperature 
guesses[1,:,:] = 3                     # Excitation  Temp
guesses[2,:,:] = 14.5                  # log(column)
guesses[3,:,:] = moment2  # Line width / 5 (the NH3 moment overestimates linewidth)               guesses[4,:,:] = moment1  # Line centroid              
guesses[5,:,:] = 0.5                   # F(ortho) - ortho NH3 fraction (fixed)

F=False
T=True


cubes.fiteach(fittype='ammonia', multifit=True, guesses=guesses,
              integral=False, verbose_level=3, fixed=[F,F,F,F,F,T], signal_cut=3,
              limitedmax=[F,F,F,F,T,T],
              maxpars=[0,0,0,0,15,1],
              limitedmin=[T,T,F,T,T,T],
              minpars=[2.73,2.73,0,0,-5,0],
              use_nearest_as_guess=True, start_from_point=(xmax,ymax),
              multicore=1,
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
fitcubefile.writeto("parameter_maps.fits",clobber=True)

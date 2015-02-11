import pyspeckit
import astropy.io.fits as fits
import numpy as np
import os
from spectral_cube import SpectralCube
import signal_id
from radio_beam import Beam

OneOneIntegrated = 'NGC1333_w11.fits'
OneOneFile = ''
TwoTwoFile = ''
ThreeThreeFile = ''

w11 = fits.getdata(OneOneIntegrated)

threshold = 0.2
mask = w11.squeeze() > threshold
cube11 = pyspeckit.Cube(OneOneFile,maskmap=mask)
cube11.units="K"
cube22 = pypseckit.Cube(TwoTwoFile,maskmap=mask)
cube22.units="K"
cube33 = pypseckit.Cube(ThreeThreeFile,maskmap=mask)
cube33.units="K"

beam11 = Beam.from_fits_header(fits.getheader(OneOneFile))
cube11sc = SpectralCube.read(OneOneFile)
noisy = signal_id.noise(cube11sc, beam = beam11)

cubes = pyspeckit.CubeStack([cube11,cube22,cube33],maskmap=mask)
cubes.units="K"

peakloc = np.argmax(w11)
xmax,ymax = np.unravel_index(peakloc,w11.shape)

guesses = np.zeros((6,)+cubes.cube.shape[1:])
guesses[0,:,:] = 20                    # Kinetic temperature 
guesses[1,:,:] = 5                     # Excitation  Temp
guesses[2,:,:] = 14.5                  # log(column)
guesses[3,:,:] = cube11sc(2,0)         # Line width / 5 (the NH3 moment overestimates linewidth)                  
guesses[4,:,:] = cube11sc.mom(1,0)     # Line centroid              
guesses[5,:,:] = 0.5                   # F(ortho) - ortho NH3 fraction (fixed)

cubes.fiteach(fittype='ammonia', multifit=True, guesses=guesses,
              integral=False, verbose_level=3, fixed=[F,F,F,F,F,T], signal_cut=3,
              limitedmax=[F,F,F,F,T,T],
              maxpars=[0,0,0,0,101,1],
              limitedmin=[T,T,F,F,T,T],
              minpars=[2.73,2.73,0,0,91,0],
              use_nearest_as_guess=True, start_from_point=(xmax,ymax),
              multicore=1,
              errmap=errmap11)

fitcubefile = pyfits.PrimaryHDU(data=np.concatenate([cubes.parcube,cubes.errcube]), header=cubes.header)
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
fitcubefile.writeto("parameter_maps.fits")

import os
import gridregion

from astropy.io import fits

os.chdir('/lustre/pipeline/scratch/GAS/images/B18')

gridregion.griddata( rootdir='/lustre/pipeline/scratch/GAS', region='B18',dirname='B18_NH3_11')

hd_temp=fits.getheader('B18_NH3_11.fits')
gridregion.griddata( rootdir='/lustre/pipeline/scratch/GAS', region='B18',dirname='B18_NH3_22', templateHeader=hd_temp)
gridregion.griddata( rootdir='/lustre/pipeline/scratch/GAS', region='B18',dirname='B18_NH3_33', templateHeader=hd_temp)
gridregion.griddata( rootdir='/lustre/pipeline/scratch/GAS', region='B18',dirname='B18_C2S', templateHeader=hd_temp)
gridregion.griddata( rootdir='/lustre/pipeline/scratch/GAS', region='B18',dirname='B18_HC5N', templateHeader=hd_temp)
gridregion.griddata( rootdir='/lustre/pipeline/scratch/GAS', region='B18',dirname='B18_HC7N_21_20', templateHeader=hd_temp)
gridregion.griddata( rootdir='/lustre/pipeline/scratch/GAS', region='B18',dirname='B18_HC7N_22_21', templateHeader=hd_temp)

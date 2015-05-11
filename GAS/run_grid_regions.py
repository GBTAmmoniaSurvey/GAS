
from . import gridregion

#from astropy.io import fits


def grid_SerAqu():
	gridregion.griddata( rootdir='/lustre/pipeline/scratch/GAS', region='Serpens_Aquila',dirname='Serpens_Aquila_NH3_11')

	hd_temp=fits.getheader('Serpens_Aquila_NH3_11.fits')
	gridregion.griddata( rootdir='/lustre/pipeline/scratch/GAS', region='Serpens_Aquila',dirname='Serpens_Aquila_NH3_22', templateHeader=hd_temp)
	gridregion.griddata( rootdir='/lustre/pipeline/scratch/GAS', region='Serpens_Aquila',dirname='Serpens_Aquila_NH3_33', templateHeader=hd_temp)
	gridregion.griddata( rootdir='/lustre/pipeline/scratch/GAS', region='Serpens_Aquila',dirname='Serpens_Aquila_C2S', templateHeader=hd_temp)
	gridregion.griddata( rootdir='/lustre/pipeline/scratch/GAS', region='Serpens_Aquila',dirname='Serpens_Aquila_HC5N', templateHeader=hd_temp)
	gridregion.griddata( rootdir='/lustre/pipeline/scratch/GAS', region='Serpens_Aquila',dirname='Serpens_Aquila_HC7N_21_20', templateHeader=hd_temp)
	gridregion.griddata( rootdir='/lustre/pipeline/scratch/GAS', region='Serpens_Aquila',dirname='Serpens_Aquila_HC7N_22_21', templateHeader=hd_temp)


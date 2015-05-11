
from . import gridregion

#from astropy.io import fits

def grid_SerAqu():
    print("You will image the GBT Ammonia Survey data for Serpens_Aquila")
    data_dir='/lustre/pipeline/scratch/GAS'
    region_name='Serpens_Aquila'
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_11')
    
    hd_temp=fits.getheader(region_name'_NH3_11.fits')
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_22',     templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_33',     templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_C2S',        templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC5N',       templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_21_20', templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_22_21', templateHeader=hd_temp)

def grid_NGC1333():
    print("You will image the GBT Ammonia Survey data for NGC1333")
    data_dir='/lustre/pipeline/scratch/GAS'
    region_name='NGC1333'
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_11')
    
    hd_temp=fits.getheader(region_name'_NH3_11.fits')
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_22',     templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_33',     templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_C2S',        templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC5N',       templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_21_20', templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_22_21', templateHeader=hd_temp)

def grid_L1455():
    print("You will image the GBT Ammonia Survey data for L1455")
    data_dir='/lustre/pipeline/scratch/GAS'
    region_name='L1455'
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_11')
    
    hd_temp=fits.getheader(region_name'_NH3_11.fits')
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_22',     templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_33',     templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_C2S',        templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC5N',       templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_21_20', templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_22_21', templateHeader=hd_temp)

def grid_L1688():
    print("You will image the GBT Ammonia Survey data for L1688")
    data_dir='/lustre/pipeline/scratch/GAS'
    region_name='L1688'
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_11')
    
    hd_temp=fits.getheader(region_name'_NH3_11.fits')
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_22',     templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_33',     templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_C2S',        templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC5N',       templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_21_20', templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_22_21', templateHeader=hd_temp)

def grid_OrionA():
    print("You will image the GBT Ammonia Survey data for OrionA")
    data_dir='/lustre/pipeline/scratch/GAS'
    region_name='OrionA'
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_11')
    
    hd_temp=fits.getheader(region_name'_NH3_11.fits')
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_22',     templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_33',     templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_C2S',        templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC5N',       templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_21_20', templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_22_21', templateHeader=hd_temp)

def grid_OrionB():
    print("You will image the GBT Ammonia Survey data for OrionB")
    data_dir='/lustre/pipeline/scratch/GAS'
    region_name='OrionB'
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_11')
    
    hd_temp=fits.getheader(region_name'_NH3_11.fits')
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_22',     templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_33',     templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_C2S',        templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC5N',       templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_21_20', templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_22_21', templateHeader=hd_temp)

def grid_B18():
    print("You will image the GBT Ammonia Survey data for B18")
    data_dir='/lustre/pipeline/scratch/GAS'
    region_name='B18'
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_11')
    
    hd_temp=fits.getheader(region_name'_NH3_11.fits')
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_22',     templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_33',     templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_C2S',        templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC5N',       templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_21_20', templateHeader=hd_temp)
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_22_21', templateHeader=hd_temp)


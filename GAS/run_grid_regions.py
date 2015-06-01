
from . import gridregion

from astropy.io import fits

def grid_SerAqu():
    print("You will image the GBT Ammonia Survey data for Serpens_Aquila")
    data_dir = '/lustre/pipeline/scratch/GAS'
    region_name = 'Serpens_Aquila'
    startChannel = 1024 + 668 # default 1024
    endChannel = 1024 + 1452  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel)
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')
    startChannel = 1024 + 630 # default 1024
    endChannel = 1024 + 1452  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_22', 
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 790 # default 1024
    endChannel = 1024 + 1290  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_33', 
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 790 # default 1024
    endChannel = 1024 + 1290  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_C2S', 
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 790 # default 1024
    endChannel = 1024 + 1290  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC5N', 
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 950 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_21_20', 
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 950 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_22_21', 
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)

def grid_NGC1333():
    print("You will image the GBT Ammonia Survey data for NGC1333")
    data_dir='/lustre/pipeline/scratch/GAS'
    region_name='NGC1333'
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel)
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')
    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_C2S', 
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC5N', 
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 880
    endChannel = 1024 + 1160
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_21_20', 
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_22_21', 
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)

def grid_L1455():
    print("You will image the GBT Ammonia Survey data for L1455")
    data_dir='/lustre/pipeline/scratch/GAS'
    region_name='L1455'
    startChannel = 1024 + 630
    endChannel = 1024 + 1380
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel)
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')
    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 700 # No lines. Using the same as NH3_22 
    endChannel = 1024 + 1350 # No lines. Using the same as NH3_22
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 700 
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 700 # No Lines. Using the same as previous
    endChannel = 1024 + 1350 # No Lines. Using the same as previous
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 700 + 180 # No Lines. Reduce channel range to avoid absorption from frequency switching
    endChannel   = 1024 + 700 + 460 # No Lines. Reduce channel range to avoid absorption from frequency switching
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 700 # No Lines. Using the same as previous
    endChannel = 1024 + 1350 # No Lines. Using the same as previous
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)

def grid_L1688():
    print("You will image the GBT Ammonia Survey data for L1688")
    data_dir='/lustre/pipeline/scratch/GAS'
    region_name='L1688'
    startChannel = 1024 + 0 #As is in the fits file
    endChannel = 1024 + 650
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_11',
        startChannel = startChannel, endChannel = endChannel )
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')
    startChannel = 1024 + 0 #As is in the fits file
    endChannel = 1024 + 650
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 0 #As is in the fits file
    endChannel = 1024 + 650
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 0 #As is in the fits file
    endChannel = 1024 + 650
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 0 #As is in the fits file
    endChannel = 1024 + 650
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 0 #As is in the fits file
    endChannel = 1024 + 650
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 0 #As is in the fits file
    endChannel = 1024 + 650
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)

def grid_OrionA():
    print("You will image the GBT Ammonia Survey data for OrionA")
    data_dir='/lustre/pipeline/scratch/GAS'
    region_name='OrionA'
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_11',
        startChannel = startChannel, endChannel = endChannel)
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')
    startChannel = 1024 + 596 # default 1024
    endChannel = 1024 + 1470  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 577 # default 1024
    endChannel = 1024 + 1540  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 790 # default 1024
    endChannel = 1024 + 1290  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC5N', 
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1130  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_21_20', 
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_22_21', 
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)

def grid_OrionB():
    print("You will image the GBT Ammonia Survey data for OrionB")
    data_dir='/lustre/pipeline/scratch/GAS'
    region_name='OrionB'
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_11')
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')
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
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_11',
        startChannel = startChannel, endChannel = endChannel)
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')
    startChannel = 1024 + 596 # default 1024
    endChannel = 1024 + 1470  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 577 # default 1024
    endChannel = 1024 + 1540  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 790 # default 1024
    endChannel = 1024 + 1290  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)

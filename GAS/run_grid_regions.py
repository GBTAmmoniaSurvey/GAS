
from . import gridregion
import textwrap
import sys
import os

from astropy.io import fits

quit_message=textwrap.dedent("""\
    Release parameters not defined. This region is either not
    processed in this release or it is not yet implemented.""")

info_message="You will image the GBT Ammonia Survey data for "

gbt_dir='/lustre/pipeline/scratch/GAS'
#dr1_dir=os.getcwd()+'/lustre/pipeline/scratch/GAS'
dr1_dir='/lustre/pipeline/scratch/GAS'

def grid_SerAqu(release=None):
    """
    Function to image the Serpens_Aquila data. The release parameter is used 
    to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR1':
        file_extension='_DR1'
        mySessions = None
        data_dir=dr1_dir
    else:
        sys.exit(quit_message)
    region_name = 'Serpens_Aquila'
    print(info_message+region_name)
    #
    startChannel = 1024 + 668 # default 1024
    endChannel = 1024 + 1452  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')

    startChannel = 1024 + 630 # default 1024
    endChannel = 1024 + 1452  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 790 # default 1024
    endChannel = 1024 + 1290  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 790 # default 1024
    endChannel = 1024 + 1290  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 790 # default 1024
    endChannel = 1024 + 1290  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 950 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 950 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

def grid_NGC1333(release=None):
    """
    Function to image the NGC1333 data. The release parameter is used to 
    select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR1':
        file_extension='_DR1'
        mySessions = None
        data_dir=dr1_dir
    else:
        sys.exit(quit_message)
    region_name='NGC1333'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')

    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_C2S', 
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC5N', 
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 880
    endChannel = 1024 + 1160
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_21_20', 
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_22_21', 
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

def grid_L1455(release=None):
    """
    Function to image the L1455 data. The release parameter is used to 
    select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR1':
        file_extension='_DR1'
        mySessions = None
        data_dir=dr1_dir
    else:
        sys.exit(quit_message)
    region_name = 'L1455'
    print(info_message+region_name)
    #
    startChannel = 1024 + 630
    endChannel = 1024 + 1380
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')
    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 700 # No lines. Using the same as NH3_22 
    endChannel = 1024 + 1350 # No lines. Using the same as NH3_22
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 700 
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 700 # No Lines. Using the same as previous
    endChannel = 1024 + 1350 # No Lines. Using the same as previous
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 700 + 180 # No Lines. Reduce channel range to avoid absorption from frequency switching
    endChannel   = 1024 + 700 + 460 # No Lines. Reduce channel range to avoid absorption from frequency switching
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 700 # No Lines. Using the same as previous
    endChannel = 1024 + 1350 # No Lines. Using the same as previous
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

def grid_L1451(release=None):
    """
    Function to image the L1451 data. The release parameter is used 
    to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR1':
        file_extension='_DR1'
        mySessions = None
        data_dir=dr1_dir
    else:
        sys.exit(quit_message)
    region_name = 'L1451'
    print(info_message+region_name)
    #
    startChannel = 1024 + 630
    endChannel = 1024 + 1380
        gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')
    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 700 # No lines. Using the same as NH3_22 
    endChannel = 1024 + 1350 # No lines. Using the same as NH3_22
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 700 
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 700 # No Lines. Using the same as previous
    endChannel = 1024 + 1350 # No Lines. Using the same as previous
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 700 + 180 # No Lines. Reduce channel range to avoid absorption from frequency switching
    endChannel   = 1024 + 700 + 460 # No Lines. Reduce channel range to avoid absorption from frequency switching
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 700 # No Lines. Using the same as previous
    endChannel = 1024 + 1350 # No Lines. Using the same as previous
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

def grid_L1688(release=None):
    """
    Function to image the L1688 data. The release parameter is used 
    to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR1':
        file_extension='_DR1'
        mySessions = None
        data_dir=dr1_dir
    else:
        sys.exit(quit_message)
    region_name = 'L1688'
    print(info_message+region_name)
    #
    startChannel = 1024 + 650 #As is in the fits file
    endChannel = 1024 + 1400
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)

    hd_temp=fits.getheader(region_name+'_NH3_11.fits')

    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 700
    endChannel = 1024 + 1350
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

def grid_OrionA(release=None):
    """
    Function to image the OrionA data. The release parameter is used 
    to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR1':
        file_extension='_DR1'
        mySessions = None
        data_dir=dr1_dir
    else:
        sys.exit(quit_message)
    region_name = 'OrionA'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')

    startChannel = 1024 + 596 # default 1024
    endChannel = 1024 + 1470  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 577 # default 1024
    endChannel = 1024 + 1540  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 790 # default 1024
    endChannel = 1024 + 1290  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1130  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

def grid_OrionB(release=None):
    """
    Function to image the OrionB data. The release parameter is used 
    to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR1':
        file_extension='_DR1'
        mySessions = None
        data_dir=dr1_dir
    else:
        sys.exit(quit_message)
    region_name = 'OrionB'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)

    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')

    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

def grid_B18(release=None):
    """
    Function to image the B18 data. The release parameter is used 
    to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR1':
        file_extension='_DR1'
        mySessions = None
        data_dir=dr1_dir
    else:
        sys.exit(quit_message)
    region_name = 'B18'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')

    startChannel = 1024 + 596 # default 1024
    endChannel = 1024 + 1470  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 577 # default 1024
    endChannel = 1024 + 1540  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 790 # default 1024
    endChannel = 1024 + 1290  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)
    
    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

def grid_Cepheus(release=None):
    """
    Function to image the Cepheus data. The release parameter is used 
    to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR1':
        file_extension='_DR1'
        mySessions = None
        data_dir=dr1_dir
    else:
        sys.exit(quit_message)
    region_name = 'Cepheus'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')

    startChannel = 1024 + 596 # default 1024
    endChannel = 1024 + 1470  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 577 # default 1024
    endChannel = 1024 + 1540  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 790 # default 1024
    endChannel = 1024 + 1290  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

def grid_IC5146(release=None):
    """
    Function to image the IC5146 data. The release parameter is used 
    to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR1':
        file_extension='_DR1'
        mySessions = None
        data_dir=dr1_dir
    else:
        sys.exit(quit_message)
    region_name = 'IC5146'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')

    startChannel = 1024 + 596 # default 1024
    endChannel = 1024 + 1470  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 577 # default 1024
    endChannel = 1024 + 1540  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 790 # default 1024
    endChannel = 1024 + 1290  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

def grid_B1(release=None):
    """
    Function to image the B1 data. The release parameter is used 
    to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR1':
        file_extension='_DR1'
        mySessions = None
        data_dir=dr1_dir
    else:
        sys.exit(quit_message)
    region_name = 'B1'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')

    startChannel = 1024 + 596 # default 1024
    endChannel = 1024 + 1470  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 577 # default 1024
    endChannel = 1024 + 1540  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 790 # default 1024
    endChannel = 1024 + 1290  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

def grid_IC348(release=None):
    """
    Function to image the IC348 data. The release parameter is used 
    to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR1':
        file_extension='_DR1'
        mySessions = None
        data_dir=dr1_dir
    else:
        sys.exit(quit_message)
    region_name = 'IC348'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')

    startChannel = 1024 + 596 # default 1024
    endChannel = 1024 + 1470  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 577 # default 1024
    endChannel = 1024 + 1540  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 790 # default 1024
    endChannel = 1024 + 1290  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)


def grid_B59(release=None):
    """
    Function to image the B59 data. The release parameter is used 
    to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR1':
        file_extension='_DR1'
        mySessions = None
        data_dir=dr1_dir
    else:
        sys.exit(quit_message)
    region_name = 'B59'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11.fits')

    startChannel = 1024 + 596 # default 1024
    endChannel = 1024 + 1470  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 577 # default 1024
    endChannel = 1024 + 1540  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 790 # default 1024
    endChannel = 1024 + 1290  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)

    startChannel = 1024 + 600 # default 1024
    endChannel = 1024 + 1150  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, 
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension)


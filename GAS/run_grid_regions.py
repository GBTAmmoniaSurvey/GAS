
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
dr2_dir='/lustre/pipeline/scratch/GAS'

def grid_SerAqu(release=None):
    """
    Function to image the Serpens_Aquila data. The release parameter is 
    used to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR2':
        file_extension='_DR2'
        mySessions = None
        data_dir=dr2_dir
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
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')

    startChannel = 1024 + 630 # default 1024
    endChannel = 1024 + 1452  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_22',
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

def grid_NGC1333(release=None,rbflag=True):
    """    Function to image the NGC1333 data. The release parameter is used to 
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
    # startChannel = 1024 + 655 # default 1024
    # endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', outerwindow=20, 
        Sessions=mySessions, file_extension=file_extension,
                         rebase=rbflag)
    #
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')

    # startChannel = 1724
    # endChannel = 2374
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_22', outerwindow=20,
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension,
                         rebase=rbflag)
    #
    # startChannel = 1888
    # endChannel = 2220
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_33', outerwindow=20,
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension,
                         rebase=rbflag)

	#
    # startChannel = 1724 + 165
    # endChannel = 1724 + 497
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_C2S', outerwindow=8,
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension,
                         rebase=rbflag)
    #
    # startChannel=1724 + 168
    # endChannel=1724 + 500
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC5N', 
        templateHeader=hd_temp, outerwindow=8,
        Sessions=mySessions, file_extension=file_extension,
                         rebase=rbflag)
    #
    # startChannel=1878
    # endChannel=2210
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_21_20', outerwindow=8,
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension,
                         rebase=rbflag)
    #
    # startChannel= 1890
    # endChannel= 2222
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_HC7N_22_21', outerwindow=8,
        templateHeader=hd_temp,
        Sessions=mySessions, file_extension=file_extension,
                         rebase=rbflag)

def grid_L1455(release=None):
    """
    Function to image the L1455 data. The release parameter is 
    used to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR2':
        file_extension='_DR2'
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
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')
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


def grid_L1448(release=None):
    """
    Function to image the L1448 data. The release parameter is 
    used to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR2':
        file_extension='_DR2'
        mySessions = None
        data_dir=dr1_dir
    else:
        sys.exit(quit_message)
    region_name = 'L1448'
    print(info_message+region_name)
    #
    startChannel = 1024 + 630
    endChannel = 1024 + 1380
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
    Function to image the L1451 data. The release parameter is 
    used to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR2':
        file_extension='_DR2'
        mySessions = None
        data_dir=dr2_dir
    else:
        sys.exit(quit_message)
    region_name='L1451'

    print(info_message+region_name)
    #
    startChannel = 1024 + 630
    endChannel = 1024 + 1380
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
    Function to image the L1688 data. The release parameter is 
    used to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR1':
        file_extension='_DR1'
        mySessions = range(1,45)
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
    Function to image the OrionA data. The release parameter is 
    used to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR1':
        file_extension='_DR1'
        mySessions = range(27)
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
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')

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

def grid_OrionB_NGC2023_2024(release=None):
    """
    Function to image the OrionB_NGC2068-2071 data. The release parameter is 
    used to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR2':
        file_extension='_DR2'
        mySessions = None
        data_dir=dr2_dir
    else:
        sys.exit(quit_message)
    region_name='OrionB_NGC2023-2024'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')

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

def grid_OrionB_NGC2068_2071(release=None):
    """
    Function to image the OrionB_NGC2068-2071 data. The release parameter is 
    used to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR2':
        file_extension='_DR2'
        mySessions = None
        data_dir=dr2_dir
    else:
        sys.exit(quit_message)
    region_name='OrionB_NGC2068-2071'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')
    
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
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')

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

def grid_Cepheus_L1228(release=None):
    """
    Function to image the Cepheus L1228 data. The release parameter is 
    used to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR2':
        file_extension='_DR2'
        mySessions = None
        data_dir=dr2_dir
    else:
        sys.exit(quit_message)
    region_name='Cepheus_L1228'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')

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

def grid_Cepheus_L1251(release=None):
    """
    Function to image the Cepheus L1251 data. The release parameter is 
    used to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR2':
        file_extension='_DR2'
        mySessions = None
        data_dir=dr2_dir
    else:
        sys.exit(quit_message)
    region_name='Cepheus_L1251'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')

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
    Function to image the IC5146 data. The release parameter is 
    used to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR2':
        file_extension='_DR2'
        mySessions = None
        data_dir=dr2_dir
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
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')

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
    elif release == 'DR2':
        file_extension='_DR2'
        mySessions = None
        data_dir=dr2_dir
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
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')

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
    elif release == 'DR2':
        file_extension='_DR2'
        mySessions = None
        data_dir=dr2_dir
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
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')

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
    elif release == 'DR2':
        file_extension='_DR2'
        mySessions = None
        data_dir=dr2_dir
    else:
        sys.exit(quit_message)
    region_name = 'B59'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024	# 763 output channels
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_11',
        startChannel = startChannel, endChannel = endChannel)

    hd_temp=fits.getheader(region_name+'_NH3_11.fits')
    startChannel = 1024 + 851 # default 1024	
    endChannel = 1024 + 1201  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_22',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 851 # default 1024
    endChannel = 1024 + 1201  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_NH3_33',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 847 # default 1024
    endChannel = 1024 + 1197  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_C2S',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 850 # default 1024
    endChannel = 1024 + 1200  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC5N',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 821 # default 1024
    endChannel = 1024 + 1171  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_21_20',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)
    startChannel = 1024 + 813 # default 1024
    endChannel = 1024 + 1163  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, dirname=region_name+'_HC7N_22_21',
        startChannel = startChannel, endChannel = endChannel, templateHeader=hd_temp)

def grid_HC2(release=None):
    """
    Function to image the Heiles Cloud2 data. The release parameter is used 
    to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR2':
        file_extension='_DR2'
        mySessions = None
        data_dir=dr2_dir
    else:
        sys.exit(quit_message)
    region_name = 'HC2'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')

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

def grid_L1689(release=None):
    """
    Function to image the L1689 data. The release parameter is used 
    to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR2':
        file_extension='_DR2'
        mySessions = None
        data_dir=dr2_dir
    else:
        sys.exit(quit_message)
    region_name = 'L1689'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')

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

def grid_B1E(release=None):
    """
    Function to image the B1E data. The release parameter is used 
    to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR2':
        file_extension='_DR2'
        mySessions = None
        data_dir=dr1_dir
    else:
        sys.exit(quit_message)
    region_name = 'B1E'
    print(info_message+region_name)
    #
    startChannel = 1024 + 655 # default 1024
    endChannel = 1024 + 1418  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')

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

def grid_SerMWC(release=None):
    """
    Function to image the Serpens MWC 297 data. The release parameter is 
    used to select the proper sessions to be imaged and the pre-defined 
    file extension.
    """
    if not release:
        file_extension='_all'
        mySessions=None
        data_dir=gbt_dir
    elif release == 'DR2':
        file_extension='_DR2'
        mySessions = None
        data_dir=dr2_dir
    else:
        sys.exit(quit_message)
    region_name = 'Serpens_MWC297'
    print(info_message+region_name)
    #
    startChannel = 1024 + 668 # default 1024
    endChannel = 1024 + 1452  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_11', 
        startChannel = startChannel, endChannel = endChannel, 
        Sessions=mySessions, file_extension=file_extension)
    
    hd_temp=fits.getheader(region_name+'_NH3_11'+file_extension+'.fits')

    startChannel = 1024 + 630 # default 1024
    endChannel = 1024 + 1452  # default 3072
    gridregion.griddata( rootdir=data_dir, region=region_name, 
        dirname=region_name+'_NH3_22',
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

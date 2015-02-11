#/bin/env python

import os 
import subprocess
import glob

doMap = False
SessionNumber = 1
StartScan = 11
EndScan = 58
Region = 'Perseus_map_NGC1333-A'
Window = '3'

RawDataDir = '/lustre/pipeline/scratch/GAS/rawdata/'
SessionDir = 'AGBT15A_430_'+str(SessionNumber).zfill(2)+'.raw.vegas/'
BankNames = ['A','B','C','D','E','F','G','H']

OptionDict = {'--window':Window,
              '--imaging-off':'',
              '--clobber':'',
              '-v':'4',
              '-m':'{0}:{1}'.format(StartScan,EndScan),
              '--units':'tmb',
              '--keep-temporary-files':'',
              }

# Grind through calibration
for bank in BankNames:
    InputFile = RawDataDir+SessionDir+'AGBT15A_430_'+\
        str(SessionNumber).zfill(2)+\
        '.raw.vegas.{0}.fits'.format(bank)
    command = 'gbtpipeline-test -i '+InputFile
    for key in OptionDict:
        command = command+' '+key+' '+OptionDict[key]
    print(command)
    subprocess.call(command,shell=True)

# Map to SDFITS in AIPS (we're getting monotonically farther from good)
if doMap:
    filelist = glob.glob('*scan*{0}*{1}*fits'.format(StartScan,EndScan))

    sdfList = []
    for ctr,fl in enumerate(filelist):
        sdfFile = '{0}.{1}.sdf '.format(Region,ctr)
        sdfList = sdfList + [sdfFile]
        command = 'idlToSdfits -l -o '+sdfFile+' '+fl
        print(command)
        subprocess.call(command,shell=True)

    # Contruct the database
    command = 'doImage /home/gbtpipeline/release/contrib/dbcon.py {0} '.format(os.getuid())
    for fl in sdfList:
        command = command + fl
    print(command)
    subprocess.call(command,shell=True)

    # Run the imaging
    command = 'doImage /home/gbtpipeline/release/contrib/mapDefault.py {0}'.format(os.getuid())
    print(command)
    subprocess.call(command,shell=True)

    # clear the AIPS catalog
    command = 'doImage /home/gbtpipeline/release/contrib/clear_AIPS_catalog.py {0}'.format(os.getuid())
    subprocess.call(command,shell=True)

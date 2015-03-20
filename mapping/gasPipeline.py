#/bin/env python
import pdb
import os 
import subprocess
import glob
import warnings
# doMap = False
# SessionNumber = 1
# StartScan = 11
# EndScan = 58
# Region = 'Perseus_map_NGC1333-A'
# Window = '3'

def wrapper(logfile='ObservationLog.csv',region='NGC1333',
            window=['0','1','3','4','5','6']):
    """
    This is the GAS pipeline which chomps the observation logs and
    then batch calibrates the data.  It requires AstropPy because
    their tables are pretty nifty.

    Usage from iPython prompt: In the directory you want to produce
    calibrated data in source this file execfile() and then run the
    script wrapper:

    execfile('/users/erosolow/GAS/mapping/gasPipeline.py')
    wrapper(logfile='../ObservationLog.csv',region='NGC1333',window=['3'])

    logfile -- Full path to CSV version of the logfile
    region -- Substring to calibrate a region
    window -- List of spectral windows to calibrate.
    """

    if ~os.access(logfile,os.R_OK):
        updateLogs()

    t = parseLog(logfile=logfile)
    for observation in t:
        if region == observation['Region name']:
            for thisWindow in window:
                doPipeline(SessionNumber=observation['Session'],
                           StartScan=observation['Start Scan'],
                           EndScan=observation['End Scan'],
                           Source=observation['Source'],
                           Region=region,
                           Window=str(thisWindow))


def parseLog(logfile='ObservationLog.csv'):
    try:
        from astropy.table import Table
    except:
        warnings.warn('GAS Pipeline requires astropy.  Try Anaconda!')
        return
    t = Table.read(logfile)
    return(t)

def updateLogs(output='ObservationLog.csv'):
    command = "wget --no-check-certificate --output-document="+output+" 'https://docs.google.com/spreadsheet/ccc?key=1F6MnXjK1Y1VyM8zWW3R5VvLAFF2Hkc85SGBRBxQ24JY&output=csv'"
    subprocess.call(command,shell=True)

def doPipeline(SessionNumber=1,StartScan = 11, EndScan=58, 
               Source='Perseus_map_NGC1333-A', Window='0', 
               Region = 'NGC1333', OptionDict = None,
               OutputRoot = None):

    RawDataDir = '/lustre/pipeline/scratch/GAS/rawdata/'
    SessionDir = 'AGBT15A_430_'+str(SessionNumber).zfill(2)+'.raw.vegas/'
    BankNames = ['A','B','C','D','E','F','G','H']

    WindowDict = {'0':'NH3_11',
                  '1':'HC7N_21_20',
                  '2':'C2S',
                  '3':'NH3_22',
                  '4':'NH3_33',
                  '5':'HC5N',
                  '6':'HC7N_22_21'}
    
    # Set default pipeline options as a dictionary
    if OptionDict is None:
        OptionDict = {'--window':Window,
                      '--imaging-off':'',
                      '--clobber':'',
                      '-v':'4',
                      '-m':'{0}:{1}'.format(StartScan,EndScan),
                      '--units':'tmb',
                      '--keep-temporary-files':'',
                      }
    if OutputRoot is None:
        OutputRoot = os.getcwd()+'/'
    # Try to make the output directory
    print('Region {0}'.format(Region))
    OutputDirectory = OutputRoot+Region+WindowDict[Window]
    if not os.access(OutputDirectory,os.W_OK):
        try:
            os.mkdir(OutputDirectory)
            print('Made directory {0}'.format(OutputDirectory))
        except:
            warnings.warn('Unable to make output directory '+OutputDirectory)
            raise

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
        # Pipeline outputs into master directory for region.
        for feed in ['0','1','2','3','4','5','6']:
            for pol in ['0','1']:
                indexname = Source+'_scan_{0}_{1}_window{2}_feed{3}_pol{4}.index'.\
                    format(StartScan,EndScan,Window,feed,pol) 
                try:
                    os.rename(indexname,OutputDirectory+'/'+indexname)
                except:
                    pass
                filename = Source+'_scan_{0}_{1}_window{2}_feed{3}_pol{4}.fits'.\
                    format(StartScan,EndScan,Window,feed,pol) 
                try:
                    os.rename(filename,OutputDirectory+'/'+filename)
                except:
                    pass

    # Map to SDFITS in AIPS (we're getting monotonically farther from good)
    # if doMap:
    #     filelist = glob.glob('*scan*{0}*{1}*fits'.format(StartScan,EndScan))

    #     sdfList = []
    #     for ctr,fl in enumerate(filelist):
    #         sdfFile = '{0}.{1}.sdf '.format(Source,ctr)
    #         sdfList = sdfList + [sdfFile]
    #         command = 'idlToSdfits -l -o '+sdfFile+' '+fl
    #         print(command)
    #         subprocess.call(command,shell=True)

    #     # Contruct the database
    #     command = 'doImage /home/gbtpipeline/release/contrib/dbcon.py {0} '.format(os.getuid())
    #     for fl in sdfList:
    #         command = command + fl
    #     print(command)
    #     subprocess.call(command,shell=True)

    #     # Run the imaging
    #     command = 'doImage /home/gbtpipeline/release/contrib/mapDefault.py {0}'.format(os.getuid())
    #     print(command)
    #     subprocess.call(command,shell=True)

    #     # clear the AIPS catalog
    #     command = 'doImage /home/gbtpipeline/release/contrib/clear_AIPS_catalog.py {0}'.format(os.getuid())
    #     subprocess.call(command,shell=True)

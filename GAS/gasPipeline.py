import pdb
import os 
import subprocess
import glob
import warnings
from astropy.time import Time
from . import catalogs


def move_files(region='Cepheus_L1251', session=81, 
               prefix='Cepheus_L1251_map_1_scan_26_49'):
    """
    Sometimes the pipeline fails to move the calibrated files into the proper 
    folder. 
    move_files(region='Cepheus_L1251', session='81', prefix='Cepheus_L1251_map_1_scan_26_49')

    region -- Region name. The files will be moved to folders like 
              region+line_name (eg NGC1333_C2S).

    session -- Integer with session number of the observations. This is 
               added to the original filename.

    prefix -- String with prefix of files to be searched for. 
    """
    folder=[ region+'_NH3_11',
             region+'_NH3_22',
             region+'_NH3_33',
             region+'_C2S',
             region+'_HC5N',
             region+'_HC7N_22_21',
             region+'_HC7N_21_20']
    window=['0','3','4','2','5','6','1']
    for i in range(len(folder)):
        file_list=glob.glob('{0}*window{1}*fits'.format(prefix,window[i]))
        if len(file_list) > 0:
            for file_i in file_list:
                os.rename( file_i, '{0}/{1}'.format( folder[i], 
                           file_i.replace('.fits', '_sess{0}.fits'.format(i))))

def fillAll(overwrite=False):

    RawDir = '/lustre/pipeline/scratch/GAS/rawdata/'
    try:
        os.chdir(RawDir)
    except OSError:
        warnings.warn("fillAll() must be run on GB machines with access to lustre")
        return

    catalogs.updateLogs(release=release)
    log = catalogs.parseLog()
    uniqSess = set(log['Session'].data.data)
    for session in uniqSess:
        if not overwrite:
            SessionName = 'AGBT15A_430_{0}'.format(session)
            OutputDir = SessionName+'.raw.vegas'
            if not os.access(OutputDir,os.W_OK):
                command = 'sdfits -backends=vegas AGBT15A_430_{0}'.format(session)
                subprocess.call(command,shell=True)
                groupchange = 'chgrp gas -R '+OutputDir
                subprocess.call(groupchange,shell=True)
                permissions = 'chmod g+rw -R '+OutputDir
                subprocess.call(permissions,shell=True)
        else:
            SessionName = 'AGBT15A_430_{0}'.format(session)
            OutputDir = SessionName+'.raw.vegas'
            subprocess.call('rm -rf '+OutputDir,shell=True)
            command = 'sdfits -backends=vegas AGBT15A_430_{0}'.format(session)
            subprocess.call(command,shell=True)
            groupchange = 'chgrp gas -R '+OutputDir
            subprocess.call(groupchange,shell=True)
            permissions = 'chmod g+rw -R '+OutputDir
            subprocess.call(permissions,shell=True)

def reduceAll(overwrite=False, release = 'all'):
    catalogs.updateLogs(release=release)
    catalogs.updateCatalog(release=release)
    RegionCatalog = catalogs.GenerateRegions()
    Log = catalogs.parseLog()
    uniqSrc = RegionCatalog['Region name']
    cwd = os.getcwd()
    for region in uniqSrc:
        if region != 'none':
            try:
                os.chdir(cwd+'/'+region)
            except OSError:
                os.mkdir(cwd+'/'+region)
                os.chdir(cwd+'/'+region)
            wrapper(region=region, overwrite = overwrite, 
                    release=release, obslog = Log)
            os.chdir(cwd)

def wrapper(logfile='ObservationLog.csv',region='NGC1333',
            window=['0','1','2','3','4','5','6'],
            overwrite=False,startdate = '2015-01-1',
            enddate='2020-12-31',release='all',obslog = None):
    """
    This is the GAS pipeline which chomps the observation logs and
    then batch calibrates the data.  It requires AstroPy because
    their tables are pretty nifty.
    
    wrapper(logfile='../ObservationLog.csv',region='NGC1333',window=['3'])

    region : string 
        Region name as given in logs
    window : list of strings
        List of spectral windows to calibrate
    logfile : string 
        Full path to CSV version of the logfile (optional)
    obslog : astropy.Table
        Table representing an already parsed observation log
    overwrite : bool
        If True, carries out calibration for files already present on disk.
    startdate : string 
        representation of date in format YYYY-MM-DD for beginning calibration 
    enddate : string 
        date in format YYYY-MM-DD for ending calibration 
    release : string
        name of column in the log file that is filled with boolean
        values indicating whether a given set of scans belongs to the data
        release.
    If a logfile or obslog isn't specified, logs will be retrieved from Google.
    """
    StartDate = Time(startdate)
    EndDate = Time(enddate)
    if not os.access(logfile,os.R_OK):
        catalogs.updateLogs(release=release)

    if obslog is None:
        t = catalogs.parseLog(logfile=logfile)
    else:
        t = obslog

    for observation in t:
        print(observation['Date'])
        ObsDate = Time(observation['Date'])
        if (region == observation['Region name']) & \
                (ObsDate >= StartDate) & (ObsDate <= EndDate) & \
                (observation[release] == 'TRUE'):
            for thisWindow in window:
                if str(observation['Beam Gains']) == '--':
                    Gains = '1,1,1,1,1,1,1,1,1,1,1,1,1,1'
                else:
                    Gains = observation['Beam Gains']
                if str(observation['Special RawDir']) == '--':
                    doPipeline(SessionNumber=observation['Session'],
                           StartScan=observation['Start Scan'],
                           EndScan=observation['End Scan'],
                           Source=observation['Source'],
                           Gains=Gains,
                           Region=region,
                           Window=str(thisWindow),overwrite=overwrite)
                else :
                        doPipeline(SessionNumber=observation['Session'],
                           StartScan=observation['Start Scan'],
                           EndScan=observation['End Scan'],
                           Source=observation['Source'],
                           Gains=Gains,
                           Region=region,
                           RawDataDir=observation['Special RawDir'],
                           Window=str(thisWindow),overwrite=overwrite)

<<<<<<< HEAD
def parseLog(logfile='ObservationLog.csv'):
    """
    Ingests a CSV log file into an astropy table
    """
    try:
        from astropy.table import Table
    except:
        warnings.warn('GAS Pipeline requires astropy.  Try Anaconda!')
        return
    t = Table.read(logfile)
    return(t)

def updateLogs(output='ObservationLog.csv',release=None):
    if release is None:
        command = "wget --no-check-certificate --output-document="+output+" 'https://docs.google.com/spreadsheet/ccc?key=1F6MnXjK1Y1VyM8zWW3R5VvLAFF2Hkc85SGBRBxQ24JY&output=csv'"
        # The returns from subprocess are the error codes from the OS
        # If 0 then it worked so we should return True
        return not subprocess.call(command,shell=True)
    if 'DR1' in release:
        from astropy.utils.data import get_pkg_data_filename
        filename = get_pkg_data_filename('data/ObservationLog_DR1.csv',
                                         package='GAS')
        command = 'cp '+filename+' ./ObservationLog.csv'
        return not subprocess.call(command,shell=True)
    else:
        warnings.warn('Updating logs failed for non-existent data release.')
        return False

def updateCatalog(output='RegionCatalog.csv',release=None):
    if release is None:
        command = "wget --no-check-certificate --output-document="+output+" 'https://docs.google.com/spreadsheets/d/140SUALscsm4Lco2WU3jDaREtUnf4jA9ZEBrMg4VAdKw/export?gid=1599734490&format=csv'"
        return not subprocess.call(command,shell=True)
    if 'DR1' in release:
        from astropy.utils.data import get_pkg_data_filename
        filename = get_pkg_data_filename('data/RegionCatalog_DR1.csv',
                                         package='GAS')
        command = 'cp '+filename+' ./'+output
        return not subprocess.call(command,shell=True)
    else:
        warnings.warn('Updating logs failed for non-existent data release.')
        return False

=======
>>>>>>> GBTAmmoniaSurvey/master
def doPipeline(SessionNumber=1,StartScan = 11, EndScan=58, 
               Source='Perseus_map_NGC1333-A', Window='0', 
               Region = 'NGC1333', OptionDict = None,
               RawDataDir = None, 
               Gains=None,
               OutputRoot = None, overwrite=False):

    if RawDataDir is None:
        RawDataDir = '/lustre/pipeline/scratch/GAS/rawdata/' 
    if Gains is None:
        Gains = '1,1,1,1,1,1,1,1,1,1,1,1,1,1'
    SessionDir = 'AGBT15A_430_'+str(SessionNumber).zfill(2)+'.raw.vegas/'
    BankNames = ['A','B','C','D','E','F','G','H']
    print('Reducing '+SessionDir)
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
                      '--smoothing-kernel-size':'0',
                      '--keep-temporary-files':'',
                      '--beam-scaling':Gains}
    if OutputRoot is None:
        OutputRoot = os.getcwd()+'/'
    # Try to make the output directory
    print('Region {0}'.format(Region))
    OutputDirectory = OutputRoot+Region+'_'+WindowDict[Window]
    if not os.access(OutputDirectory,os.W_OK):
        try:
            os.mkdir(OutputDirectory)
            print('Made directory {0}'.format(OutputDirectory))
        except:
            warnings.warn('Unable to make output directory '+OutputDirectory)
            raise

    for bank in BankNames:
        # Loop over each feed and polarization
        # we check if a pipeline call is necessary. 
        for feed in ['0','1','2','3','4','5','6']:
            for pol in ['0','1']:
                FilesIntact = True
                if not overwrite:
                    outputfile = Source+'_scan_{0}_{1}_window{2}_feed{3}_pol{4}_sess{5}.fits'.\
                        format(StartScan,EndScan,Window,feed,pol,SessionNumber) 
                    FilesIntact = FilesIntact and os.path.exists(OutputDirectory+'/'+outputfile)
                    if FilesIntact:
                        print('Data for Polarization {0} of Feed {1} appear on disk... skipping'.format(pol,feed))
                #
                if (not FilesIntact) or (overwrite):
                    InputFile = RawDataDir+SessionDir+'AGBT15A_430_'+\
                        str(SessionNumber).zfill(2)+\
                        '.raw.vegas.{0}.fits'.format(bank)
                    command = 'gbtpipeline-test -i '+InputFile
                    for key in OptionDict:
                        command = command+' '+key+' '+OptionDict[key]
                    command = command+' --feed '+feed+' --pol '+pol
                    print(command)
                    subprocess.call(command,shell=True)

                    indexname    = Source+'_scan_{0}_{1}_window{2}_feed{3}_pol{4}.index'.\
                        format(StartScan,EndScan,Window,feed,pol) 
                    outindexname = Source+'_scan_{0}_{1}_window{2}_feed{3}_pol{4}_sess{5}.index'.\
                        format(StartScan,EndScan,Window,feed,pol,SessionNumber) 
                    try:
                        os.rename(indexname,OutputDirectory+'/'+outindexname)
                    except:
                        pass
                    
                    filename   = Source+'_scan_{0}_{1}_window{2}_feed{3}_pol{4}.fits'.\
                        format(StartScan,EndScan,Window,feed,pol) 
                    outputfile = Source+'_scan_{0}_{1}_window{2}_feed{3}_pol{4}_sess{5}.fits'.\
                        format(StartScan,EndScan,Window,feed,pol,SessionNumber) 
                    try:
                        os.rename(filename,OutputDirectory+'/'+outputfile)
                        os.chown(OutputDirectory+'/'+outputfile,0774)
                    except:
                        pass

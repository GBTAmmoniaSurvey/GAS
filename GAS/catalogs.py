import subprocess
import glob
import os
import warnings
from astropy.table import Table, join
import numpy as np

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
    if (release is None) or ('all' in release):
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

def GenerateRegions(refresh=False,release='all'):

    if refresh:
        updateLogs(release=release)
        updateCatalog(release=release)
    if not os.access('ObservationLog.csv',os.R_OK):
        updateLogs(release = release)
    if not os.access('RegionCatalog.csv',os.R_OK):
        updateCatalog(release = release)
    obs = Table.read('ObservationLog.csv')
    cat = Table.read('RegionCatalog.csv')

# This takes out rows that are empty
# This needs to be done twice for some reason... 
    for idx, row in enumerate(cat):
        if not row['BoxName']:
            cat.remove_row(idx)

    for idx, row in enumerate(cat):
        if not row['BoxName']:
            cat.remove_row(idx)
    obs.rename_column('Source','BoxName')
    joincat = join(obs,cat,keys='BoxName')
    groupcat = joincat.group_by('Region name')
    min_values = groupcat.groups.aggregate(np.min)
    max_values = groupcat.groups.aggregate(np.max)
    mean_values = groupcat.groups.aggregate(np.mean)
    vavg = 0.5*(min_values['VLSR'] + max_values['VLSR'])
    vrange = max_values['VLSR']- min_values['VLSR']
    mean_values['VAVG'] = vavg
    mean_values['VRANGE'] = vrange

    return(mean_values)


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

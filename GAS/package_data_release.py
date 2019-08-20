import os
import glob
from . import catalogs

'''
Set of functions to tar groups of files for easier posting to 
Dataverse and/or other data repositories
Want to do:
1. Organize groups of files together:
    a) By region
    b) By line
2. Include a text file with list of files and what they are
'''

def tar_moment_files(regions=None,file_extension='all_rebase3',release='all'):
    if regions is None:
        RegionCatalog = catalogs.GenerateRegions(release=release)
    else:
        RegionCatalog = catalogs.GenerateRegions(release=release)
        keep = [idx for idx, row in enumerate(RegionCatalog) if row['Region name'] in regions]
        RegionCatalog = RegionCatalog[keep]
    outDir = 'tarFiles'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
        
    lineList = ['NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    for ThisRegion in RegionCatalog:
        region = ThisRegion['Region name']
        # First moment, rms, and cubes
        # NH3 (1,1), (2,2):
        for line in ['11','22']:
            print line
            files = ['{0}/{0}_NH3_{1}_{2}_mom0_QA_trim.fits'.format(region,line,file_extension),'{0}/{0}_NH3_{1}_{2}_mom0_sigma_QA.fits'.format(region,line,file_extension),'{0}/{0}_NH3_{1}_{2}_rms_QA_trim.fits'.format(region,line,file_extension),'{0}/{0}_NH3_{1}_{2}_trim.fits'.format(region,line,file_extension)]
            fileList = '{0}/{0}_NH3_{1}_fileList.txt'.format(region,line)
            with open(fileList,'w') as f:
                for item in files:
                    f.write("%s\n" % item)
            # Use list to tar files
            os.system('tar zcvf {4}/{0}_NH3_{1}_{2}_maps.tgz -T {3}'.format(region,line,file_extension,fileList,outDir))
        # Other lines
        for line in lineList:
            files = ['{0}/{0}_{1}_{2}_mom0_QA_trim.fits'.format(region,line,file_extension),'{0}/{0}_{1}_{2}_rms_QA_trim.fits'.format(region,line,file_extension),'{0}/{0}_{1}_{2}_trim.fits'.format(region,line,file_extension)]
            fileList = '{0}/{0}_{1}_fileList.txt'.format(region,line)
            with open(fileList,'w') as f:
                for item in files:
                    f.write("%s\n" % item)
            # Use list to tar files
            os.system('tar zcvf {4}/{0}_{1}_{2}_maps.tgz -T {3}'.format(region,line,file_extension,fileList,outDir))


def tar_fit_files(regions=None,file_extension='all_rebase3',release='all'):
    if regions is None:
        RegionCatalog = catalogs.GenerateRegions(release=release)
    else:
        RegionCatalog = catalogs.GenerateRegions(release=release)
        keep = [idx for idx, row in enumerate(RegionCatalog) if row['Region name'] in regions]
        RegionCatalog = RegionCatalog[keep]
    outDir = 'tarFiles'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    # Note that NH3 (3,3), HC7N 22-21 don't have N, tau files
    lineList1 = ['C2S','HC5N','HC7N_21_20']
    lineList2 = ['NH3_33','HC7N_22_21']
    for ThisRegion in RegionCatalog:
        region = ThisRegion['Region name']
        # NH3 fits
        files = ['{0}/{0}_parameter_maps_{1}_masked.fits'.format(region,file_extension)]
        fitResults = ['Vlsr','eVlsr','Sigma','eSigma','Tkin','eTkin','Tex','eTex','N','eN']
        for result in fitResults:
            files.append('{0}/parameterMaps/{0}_NH3_{1}_{2}_masked.fits'.format(region,result,file_extension))
        fileList = '{0}/{0}_NH3_parameterMaps_fileList.txt'.format(region)
        with open(fileList,'w') as f:
            for item in files:
                f.write("%s\n" % item)
        # Use list to tar files
        os.system('tar zcvf {3}/{0}_NH3_parameterMaps_{1}.tgz -T {2}'.format(region,file_extension,fileList,outDir))
        # Other lines
        fitResults = ['vlsr','eVlsr','sigv','eSigv','Tmb','eTmb','N','tau']
        for line in lineList1:
            files = ['{0}/{0}_{1}_{2}_param_cube_masked.fits'.format(region,line,file_extension)]
            for result in fitResults:
                files.append('{0}/parameterMaps/{0}_{1}_{2}_{3}_masked.fits'.format(region,line,file_extension,result))
            fileList = '{0}/{0}_{1}_parameterMaps_fileList.txt'.format(region,line)
            with open(fileList,'w') as f:
                for item in files:
                    f.write("%s\n" % item)
            # Use list to tar files
            os.system('tar zcvf {4}/{0}_{1}_parameterMaps_{2}.tgz -T {3}'.format(region,line,file_extension,fileList,outDir))
        # Lines without N, tau maps
        fitResults = ['vlsr','eVlsr','sigv','eSigv','Tmb','eTmb']
        for line in lineList2:
            files = ['{0}/{0}_{1}_{2}_param_cube_masked.fits'.format(region,line,file_extension)]
            for result in fitResults:
                files.append('{0}/parameterMaps/{0}_{1}_{2}_{3}_masked.fits'.format(region,line,file_extension,result))
            fileList = '{0}/{0}_{1}_parameterMaps_fileList.txt'.format(region,line)
            with open(fileList,'w') as f:
                for item in files:
                    f.write("%s\n" % item)
            # Use list to tar files
            os.system('tar zcvf {4}/{0}_{1}_parameterMaps_{2}.tgz -T {3}'.format(region,line,file_extension,fileList,outDir))

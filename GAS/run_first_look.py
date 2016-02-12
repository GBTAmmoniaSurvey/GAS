import os
import textwrap
import warnings
import numpy as np
from astropy.table import Table, join
import astropy.units as u
from spectral_cube import SpectralCube
from pyspeckit.spectrum.models.ammonia_constants import voff_lines_dict
from . import first_look
from . import gasPipeline

quit_message=textwrap.dedent("""\
    Release parameters not defined. This region is either not
    processed in this release or it is not yet implemented.""")

def GenerateRegions(refresh=False):

    if refresh:
        gasPipeline.updateLogs()
        gasPipeline.updateCatalog()

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
            
def FirstLook(regions=None, file_extension='_all'):

    RegionCatalog = GenerateRegions()
    if regions is None:
        RegionCatalog = GenerateRegions()
    else:
        RegionCatalog = GenerateRegions()
        keep = [idx for idx, row in enumerate(RegionCatalog) if row['Region name'] in regions]
        RegionCatalog = RegionCatalog[keep]

    for ThisRegion in RegionCatalog:
        region_name=ThisRegion['Region name']
        print("Now NH3(1,1)")

        vsys = ThisRegion['VAVG']*u.km/u.s
        throw = 2*u.km/u.s + ThisRegion['VRANGE']*u.km/u.s/2

        file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
        file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')

        voff11 = voff_lines_dict['oneone']
        try:
            s = SpectralCube.read(file_in)
            s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
            mask = np.ones(s.shape[0],dtype=np.bool)
            for deltav in voff11:
                mask*=(np.abs(s.spectral_axis-deltav*u.km/u.s) > throw)
            a_rms = (np.where(mask != np.roll(mask,1)))[0]
            b_rms = (np.where(mask != np.roll(mask,-1)))[0]
            index_rms=first_look.create_index(a_rms, b_rms)
            index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                                   s.closest_spectral_channel(vsys-3*u.km/u.s))
            first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
            first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)

        except IOError:
            warnings.warn("File not found: {0}".format(file_in))

        linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']

        for line in linelist:
            file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
            try:
                s = SpectralCube.read(file_in)
                s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
                a_rms = [s.closest_spectral_channel(vsys+2*throw),
                         s.closest_spectral_channel(vsys-throw)]
                b_rms = [s.closest_spectral_channel(vsys+throw),
                         s.closest_spectral_channel(vsys-2*throw)]
                index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                                       s.closest_spectral_channel(vsys-3*u.km/u.s))
                index_rms=first_look.create_index( a_rms, b_rms)

                file_out=file_in.replace(file_extension+'.fits',
                                         '_base'+file_extension+'.fits')
                first_look.baseline( file_in, file_out, 
                                     index_clean=index_rms, polyorder=1)
                first_look.peak_rms( file_out, index_rms=index_rms, 
                                     index_peak=index_peak)
            except IOError:
                warnings.warn("File not found {0}".format(file_in))

def FirstLook_OrionA(file_extension='_all'):
    """
    Function to create First Look products for OrionA. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name='OrionA'
    print("Now NH3(1,1)")
    a_rms = [  0, 158, 315, 428, 530, 693]
    b_rms = [ 60, 230, 327, 438, 604, 735]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(326,470)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    #
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    ## 2nd order polynomial
    # file_out=file_in.replace('.fits','_base2.fits')
    #file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    #first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now NH3(2,2)")
    a_rms = [  0, 260, 520, 730]
    b_rms = [150, 380, 610, 850]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(380,520)
    line='NH3_22'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    ## 2nd order polynomial
    #file_out=file_in.replace('.fits','_base2.fits')
    #file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    #first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now NH3(3,3)")
    a_rms = [ 10, 250, 530]
    b_rms = [210, 310, 930]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(410,540)
    line='NH3_33'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    ##2nd order polynomial
    #file_out=file_in.replace('.fits','_base2.fits')
    #file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    #first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now CCS")
    a_rms = [  0, 260]
    b_rms = [200, 490]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(220,250)
    line='C2S'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC5N")
    # HC5N channel range must be updated
    a_rms = [  0, 500]
    b_rms = [380, 545]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(400,480)
    line='HC5N'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC7N 21-20")
    # HC7N channel range must be updated
    a_rms = [  0, 160, 480]
    b_rms = [115, 360, 525]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(400,460)
    line='HC7N_21_20'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC7N 22-21")
    # HC7N channel range must be updated
    a_rms = [  0, 480]
    b_rms = [360, 525]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(400,460)
    line='HC7N_22_21'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)

def FirstLook_B18(file_extension='_all'):
    """
    Function to create First Look products for B18. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name='B18'
    print("Now NH3(1,1)")
    a_rms = [  0, 115, 280, 385, 490, 655]
    b_rms = [ 80, 230, 345, 455, 625, 760]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(352,381)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #

    print("Now NH3(2,2)")
    a_rms = [   0, 440]
    b_rms = [ 409, 870]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(420,435)
    line='NH3_22'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now NH3(3,3)")
    a_rms = [   0, 530]
    b_rms = [ 409, 960]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(410,485)
    line='NH3_33'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now CCS")
    a_rms = [   0, 245]
    b_rms = [ 210, 490]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(225,243)
    line='C2S'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC5N")
    a_rms = [  10, 435]
    b_rms = [ 409, 540]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(414,430)
    line='HC5N'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC7N_21_20")
    a_rms = [  10, 435]
    b_rms = [ 409, 540]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(412,430)
    line='HC7N_21_20'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC7N_22_21")
    a_rms = [  10, 435]
    b_rms = [ 409, 540]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(412,430)
    line='HC7N_22_21'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)


def FirstLook_L1688(file_extension='_all'):
    """
    Function to create First Look products for L1688. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name='L1688'
    print("Now NH3(1,1)")
    a_rms = [  0, 121, 290, 404, 505, 665]
    b_rms = [ 74, 239, 332, 447, 611, 749]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(350,377)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now NH3(2,2)")
    a_rms = [   0, 349]
    b_rms = [ 285, 649]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(298,342)
    line='NH3_22'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now NH3(3,3)")
    a_rms = [   0, 395]
    b_rms = [ 272, 649]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(298,342)
    line='NH3_33'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now CCS")
    a_rms = [   0, 369]
    b_rms = [ 278, 649]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(307,325)
    line='C2S'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC5N")
    a_rms = [   0, 358]
    b_rms = [ 288, 649]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(306,317)
    line='HC5N'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    #HC7N (21-20) shows an absorption feature at ~ 91 km/s (at 23.6951 GHz)
    #from its rest frequency (used 23.6879 GHz). There's no emission line.
    #Below are the channel indeces for the absorption feature.
    #a_rms = [  0, 520]
    #b_rms = [480, 650]
    #index_peak = np.arange(485,510)
    #
    #The code didn't produce the fits file for HC7N (22-21).

def FirstLook_L1689(file_extension='_all'):
    """
    Function to create First Look products for L1689. The 
    file_extension parameter is used to select the proper files to be 
    processed. 
    """
    region_name='L1689'
    print("Now NH3(1,1)")
    a_rms = [  0, 150, 310, 420, 530, 690]
    b_rms = [ 60, 230, 330, 440, 610, 760]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(340,420)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    # 
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 3.9*u.km/u.s
    throw = 5*u.km/u.s
    for line in linelist:
        file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+2*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-2*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
        first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_out, index_rms=index_rms, 
                             index_peak=index_peak)

def FirstLook_SerAqu(file_extension='_all'):
    """
    Function to create First Look products for Serpens_Aquila. The 
    file_extension parameter is used to select the proper files to be 
    processed. 
    """
    region_name='Serpens_Aquila'
    print("Now NH3(1,1)")
    a_rms = [  0, 150, 310, 420, 530, 690]
    b_rms = [ 60, 230, 330, 440, 610, 780]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(340,420)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    # 
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 6.35*u.km/u.s
    throw = 8*u.km/u.s
    for line in linelist:
        file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+2*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-2*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
        first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_out, index_rms=index_rms, 
                             index_peak=index_peak)

def FirstLook_L1455(file_extension='_all'):
    """
    Function to create First Look products for L1455. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name='L1455'
    print("Now NH3(1,1)")
    a_rms = [   0, 140, 300, 410, 520, 680]
    b_rms = [ 105, 270, 370, 480, 630, 745]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(350,430)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)

    print("Now NH3(2,2)")
    a_rms = [   0, 340]
    b_rms = [ 290, 648]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(260,400)
    line='NH3_22'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)

    print("Now NH3(3,3)")
    a_rms = [   0, 340]  # No lines. Using the same as NH3(2,2)
    b_rms = [ 290, 648]  # No lines. Using the same as NH3(2,2)
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(260,400)  # No lines. Using the same as NH3(2,2)
    line='NH3_33'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)

    print("Now CCS")
    a_rms = [   0, 350]  
    b_rms = [ 290, 648]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(309,334)
    line='C2S'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)

    print("Now HC5N")
    a_rms = [   0, 350]  
    b_rms = [ 290, 648]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(315,325)
    line='HC5N'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)

    print("Now HC7N_21_20")
    a_rms = [   0, 180]  
    b_rms = [ 130, 275]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(128,147)
    line='HC7N_21_20'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)

    print("Now HC7N_22_21")
    a_rms = [   0, 340]  # No lines. Using the same as HC7N_21_20
    b_rms = [ 290, 648]  # No lines. Using the same as HC7N_21_20
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(308,328)  # No lines. Using the same as HC7N_21_20
    line='HC7N_22_21'
    file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)


def FirstLook_NGC1333(file_extension='_all'):
    """
    Function to create First Look products for NGC1333. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name='NGC1333'
    print("Now NH3(1,1)")
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    a_rms = [  0, 158, 315, 428, 530, 693, 751]
    b_rms = [ 60, 230, 327, 438, 604, 735, 760]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(326,430)
    first_look.baseline( file_in, file_out, index_clean=index_rms, 
                                  polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)

    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 7.9*u.km/u.s
    throw = 2.0*u.km/u.s
    for line in linelist:
        file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+20*throw),
                 s.closest_spectral_channel(vsys-20*throw)]
        b_rms = [s.closest_spectral_channel(vsys+2*throw),
                 s.closest_spectral_channel(vsys-2*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+throw),
                              s.closest_spectral_channel(vsys-throw))
        index_rms=first_look.create_index( a_rms, b_rms)
        #
        file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
        first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_out, index_rms=index_rms, 
                             index_peak=index_peak)

def FirstLook_B1(file_extension='_all'):
    """
    Function to create First Look products for B1. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name='B1'
    print("Now NH3(1,1)")
    a_rms = [  0, 130, 290, 400, 500, 660]
    b_rms = [ 70, 240, 340, 440, 620, 740]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(340,400)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, 
                                  polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)

    print("Now the rest")
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 6.6*u.km/u.s
    throw = 2.0*u.km/u.s
    for line in linelist:
        file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+3*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-3*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
        first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_out, index_rms=index_rms, 
                             index_peak=index_peak)
        
def FirstLook_IC348(file_extension='_all'):
    """
    Function to create First Look products for IC348. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name='IC348'
    print("Now NH3(1,1)")
    a_rms = [  0, 130, 290, 400, 500, 660]
    b_rms = [ 70, 240, 340, 440, 620, 740]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(340,400)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, 
                                  polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    print("Now NH3(2,2)")
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 9.0*u.km/u.s
    throw = 2.0*u.km/u.s
    for line in linelist:
        file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+3*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-3*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
        first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_out, index_rms=index_rms, 
                             index_peak=index_peak)

def FirstLook_B59(file_extension='_all'):
    """
    Function to create First Look products for B59. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name='B59'
    print("Now NH3(1,1)")
    a_rms = [  0, 130, 290, 400, 500, 660]
    b_rms = [ 70, 240, 340, 440, 620, 740]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(340,400)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    print("Now NH3(2,2)")
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 3.5*u.km/u.s
    throw = 2.0*u.km/u.s
    for line in linelist:
        file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+3*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-3*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
        first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_out, index_rms=index_rms, 
                             index_peak=index_peak)
        
def FirstLook_Cepheus(file_extension='_all'):
    """
    Function to create First Look products for Cepheus. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name = 'Cepheus'
    print("Now NH3(1,1)")
    a_rms = [  0, 135, 290, 405, 505, 665]
    b_rms = [ 70, 245, 350, 455, 625, 740]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(350,410)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    print("Now NH3(2,2)")
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = -3.8*u.km/u.s
    throw = 2.0*u.km/u.s
    for line in linelist:
        file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+3*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-3*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
        first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_out, index_rms=index_rms, 
                             index_peak=index_peak)

def FirstLook_B1E(file_extension='_all'):
    """
    Function to create First Look products for B1E. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name = 'B1E'
    print("Now NH3(1,1)")
    a_rms = [  0, 135, 290, 405, 505, 665]
    b_rms = [ 70, 245, 350, 455, 625, 740]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(350,410)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    print("Now NH3(2,2)")
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 7.3*u.km/u.s
    throw = 2.0*u.km/u.s
    for line in linelist:
        file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+3*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-3*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
        first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_out, index_rms=index_rms, 
                             index_peak=index_peak)

def FirstLook_HC2(file_extension='_all'):
    """
    Function to create First Look products for Heiles cloud2. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name = 'HC2'
    print("Now NH3(1,1)")
    a_rms = [  0, 135, 290, 405, 505, 665]
    b_rms = [ 70, 245, 350, 455, 625, 740]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(350,410)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    #
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    print("Now NH3(2,2)")
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 5.3*u.km/u.s
    throw = 2.0*u.km/u.s
    for line in linelist:
        file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+3*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-3*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
        first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_out, index_rms=index_rms, 
                             index_peak=index_peak)

def FirstLook_TMC1(file_extension='_all'):
    """
    Function to create First Look products for TMC1. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name = 'TMC1'
    print("Now NH3(1,1)")
    a_rms = [  0, 135, 290, 405, 505, 665]
    b_rms = [ 70, 245, 350, 455, 625, 740]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(350,410)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    #
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    print("Now NH3(2,2)")
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 6.0*u.km/u.s
    throw = 2.0*u.km/u.s
    for line in linelist:
        file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+3*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-3*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
        first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_out, index_rms=index_rms, 
                             index_peak=index_peak)

def FirstLook_OrionB_NGC2023_2024(file_extension='_all'):
    """
    Function to create First Look products for OrionB NGC2023-2024. The 
    file_extension parameter is used to select the proper files to be 
    processed.
    """
    region_name = 'OrionB_NGC2023-2024'
    print("Now NH3(1,1)")
    a_rms = [  0, 150, 310, 420, 520, 680]
    b_rms = [ 70, 225, 325, 435, 600, 740]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(350,410)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    #
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    print("Now NH3(2,2)")
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 10.2*u.km/u.s
    throw = 2.0*u.km/u.s
    for line in linelist:
        file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+3*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-3*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
        first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_out, index_rms=index_rms, 
                             index_peak=index_peak)

def FirstLook_OrionB_NGC2068_2071(file_extension='_all'):
    """
    Function to create First Look products for OrionB_NGC2068_2071. The 
    file_extension parameter is used to select the proper files to be 
    processed.
    """
    region_name = 'OrionB_NGC2068-2071'
    print("Now NH3(1,1)")
    a_rms = [  0, 120, 270, 390, 480, 640]
    b_rms = [ 60, 230, 330, 440, 600, 740]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(330,390)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    #
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    print("Now NH3(2,2)")
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 10.0*u.km/u.s
    throw = 2.0*u.km/u.s
    for line in linelist:
        file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+3*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-3*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
        first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_out, index_rms=index_rms, 
                             index_peak=index_peak)

def FirstLook_L1451(file_extension='_all'):
    """
    Function to create First Look products for L1451. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name = 'L1451'
    print("Now NH3(1,1)")
    a_rms = [  0, 155, 310, 420, 525, 680]
    b_rms = [ 70, 245, 350, 460, 625, 740]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(350,415)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    #
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    print("Now NH3(2,2)")
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 4.3*u.km/u.s
    throw = 2.0*u.km/u.s
    for line in linelist:
        file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+3*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-3*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
        first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_out, index_rms=index_rms, 
                             index_peak=index_peak)

def FirstLook_template(file_extension='_all'):
    """
    Function to create First Look products for TEMPLATE. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name = 'TEMPLATE'
    print("Now NH3(1,1)")
    a_rms = [  0, 135, 290, 405, 505, 665]
    b_rms = [ 70, 245, 350, 455, 625, 740]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(350,410)
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    file_out=file_in.replace(file_extension+'.fits',
                             '_base'+file_extension+'.fits')
    #
    first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    print("Now NH3(2,2)")
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 7.3*u.km/u.s
    throw = 2.0*u.km/u.s
    for line in linelist:
        file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+3*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-3*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace(file_extension+'.fits',
                                 '_base'+file_extension+'.fits')
        first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_out, index_rms=index_rms, 
                             index_peak=index_peak)

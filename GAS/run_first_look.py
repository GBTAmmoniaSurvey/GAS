import os
import textwrap
import warnings
import numpy as np
from astropy.table import Table, join
import astropy.units as u
from spectral_cube import SpectralCube
from pyspeckit.spectrum.models.ammonia_constants import voff_lines_dict
from . import first_look
from . import gasPipeline
from . import catalogs
from . import baseline

def FirstLook(regions=None, file_extension=None, release='all'):
    """
    This runs through cubes in a directory tree and generates first
    look products for each of them.  This assumes a directory naming
    convention specified in our observing logs.
    ----------
    regions : list
        List of region names (strings) to be included.  If empty, all
        regions in the log file are searched for and reduced.
    release : string
        Name of data release.  Must match boolean column name in the 
        Observation Log.
    file_extension : string
        Name of file extensions to be searched for.  Defaults to release name.  

    Note: The GAS file naming convention is
    REGION_LINENAME_EXTENSION.fits.  For example, for NGC1333 in
    ammonia (1,1), this would look for 

    NGC1333_NH3_11_all.fits
    """

    if file_extension is None:
        file_extension = '_'+release

    RegionCatalog = catalogs.GenerateRegions(release=release)
    if regions is None:
        RegionCatalog = catalogs.GenerateRegions(release=release)
    else:
        RegionCatalog = catalogs.GenerateRegions(release=release)
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

def plot_all_moments(file_extension='base_all'):
    # Get list of regions - run from images/ directory
    # Assume directories correspond to regions to be imaged
    region_list = glob("*/")
    for i in range(len(region_list)):
        region_list[i] = region_list[i].strip("/")

    line_list  = ['NH3_11','NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    label_list = ['NH$_3$(1,1)','NH$_3$(2,2)','NH$_3$(3,3)','C$_2$S','HC$_5$N',
                  'HC$_7$N (21-20)','HC$_7$N (22-21)']
    extension = file_extension
    color_table='magma'
    text_color='black'
    text_size = 14
    beam_color='#d95f02'  # previously used '#E31A1C'
    # Try single set of contours for first look images
    w11_step = 0.3
    cont_levs=2**np.arange( 0,20)*w11_step

    # Masking of small (noisy) regions
    selem = np.array([[0,1,0],[1,1,1],[0,1,0]])

    for region in region_list:
        file_w11='{0}/{0}_NH3_11_{1}_mom0.fits'.format(region,extension)
        if os.path.isfile(file_w11):
            LowestContour= cont_levs[0]*0.5
            w11_hdu = fits.open(file_w11)
            map = w11_hdu[0].data
            mask = binary_opening(map > LowestContour, selem)
            MaskedMap = mask*map
            w11_hdu[0].data = MaskedMap
            for i in range(len(line_list)):
                line_i=line_list[i]
                label_i=label_list[i]
                file_mom0='{0}/{0}_{1}_{2}_mom0.fits'.format(region,line_i,extension)
                if os.path.isfile(file_mom0):
                    line_hdu = fits.open(file_mom0)
                    # Use percentiles to set initial plot colourscale ranges
                    v_min=np.nanpercentile(line_hdu[0].data,0.1)
                    v_max=np.nanpercentile(line_hdu[0].data,99.9)
                    fig=aplpy.FITSFigure(file_mom0, hdu=0)
                    if line_i == 'NH3_11':
                        fig.show_colorscale(cmap=color_table,vmin=v_min, vmax=v_max, stretch='log',
                                            vmid=v_min-(1.*np.abs(v_min)))
                                            cbar_ticks = [0,3,6,12,24,48,96]
                        # add colorbar
                        fig.add_colorbar()
                        #fig.colorbar.set_width(0.15)
                        fig.colorbar.show(box_orientation='horizontal', width=0.1, pad=0.0, ticks=cbar_ticks,
                                                          location='top', axis_label_text='Integrated Intensity (K km s$^{-1}$)')
                    elif (line_i in ['NH3_22','NH3_33']) :
                        fig.show_colorscale(cmap=color_table,vmin=v_min, vmax=v_max, stretch='linear',
                                            vmid=v_min-(1.*np.abs(v_min)),cbar_ticks = [0,1,2,3,6,12])
                        # add colorbar
                        fig.add_colorbar()
                        #fig.colorbar.set_width(0.15)
                        fig.colorbar.show(box_orientation='horizontal', width=0.1, pad=0.0, ticks= cbar_ticks,
                                          location='top', axis_label_text='Integrated Intensity (K km s$^{-1}$)')
                    else:
                        fig.show_colorscale( cmap=color_table,vmin=v_min, vmax=v_max)
                        # add colorbar
                        fig.add_colorbar()
                        #fig.colorbar.set_width(0.15)
                        fig.colorbar.show( box_orientation='horizontal', width=0.1, pad=0.0,
                                          location='top', axis_label_text='Integrated Intensity (K km s$^{-1}$)')
                    fig.colorbar.set_font(family='sans_serif',size=text_size)
                    fig.colorbar.set_axis_label_font(family='sans_serif',size=text_size)
                    fig.set_nan_color('0.95')
                    #
                    fig.show_contour(w11_hdu, colors='gray', levels=cont_levs)
                    # Axis labels
                    fig.axis_labels.set_font(family='sans_serif',size=text_size)
                    # Ticks
                    fig.ticks.set_color(text_color)
                    fig.tick_labels.set_font(family='sans_serif',size=text_size)
                    fig.tick_labels.set_style('colons')
                    fig.tick_labels.set_xformat('hh:mm:ss')
                    fig.tick_labels.set_yformat('dd:mm')
                    # Add beam
                    fig.add_beam(major=0.0088441,minor=0.0088441,angle=0)
                    fig.beam.set_color(beam_color)
                    fig.beam.set_corner('bottom left')
                    '''
                    # Scale bar
                    # magic line of code to obtain scale in arcsec obtained from
                    # http://www.astropy.org/astropy-tutorials/Quantities.html
                    ang_sep =  (plot_param['scalebar_size'].to(u.au)/plot_param['distance']).to(u.arcsec, equivalencies dimensionless_angles())
                    fig.add_scalebar(ang_sep.to(u.degree))
                    fig.scalebar.set_corner(plot_param['scalebar_pos'])
                    fig.scalebar.set_font(family='sans_serif',size=text_size)
                    fig.scalebar.set(color=text_color)
                    fig.scalebar.set_label('{0:4.2f}'.format(plot_param['scalebar_size']))
                    '''
                    # Labels
                    fig.add_label(0.025, 0.1,
                                  '{0}\n{1}'.format(region,label_i),
                                  relative=True, color=text_color,
                                  horizontalalignment='left',
                                  family='sans_serif',size=text_size)
                    # fig.set_system_latex(True)
                    fig.save( 'figures/{0}_{1}_{2}_mom0_map.pdf'.format(region,line_i,extension),adjust_bbox=True)
                    fig.close()
            else:
                print('File {0} not found'.format(file_mom0))
    else:
        print('File {0} not found'.format(file_w11))


def FirstLook_OrionA(file_extension='_all'):
    """
    Function to create First Look products for OrionA. The file_extension 
    parameter is used to select the proper files to be processed. 
    """

    region_name='OrionA'
    vsys = 10.*u.km/u.s
    throw = 4.0*u.km/u.s

    print("Now NH3(1,1)")
    
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)

    s = SpectralCube.read(file_in)
    s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
    spaxis = s.spectral_axis.value
    index_rms = baseline.ammoniaWindow(spaxis,spaxis,window=4,v0=vsys.value)
    index_peak= ~baseline.tightWindow(spaxis,spaxis,window=3,v0=vsys.value)
    first_look.peak_rms( file_in, index_rms=index_rms, index_peak=index_peak)

    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']

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
        first_look.peak_rms( file_in, index_rms=index_rms, 
                             index_peak=index_peak)




    #region_name='OrionA'
    #print("Now NH3(1,1)")
    #a_rms = [  0, 158, 315, 428, 530, 693]
    #b_rms = [ 60, 230, 327, 438, 604, 735]
    #index_rms=first_look.create_index( a_rms, b_rms)
    #index_peak=np.arange(326,470)
    #file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    #file_out=file_in.replace(file_extension+'.fits',
    #                         '_base'+file_extension+'.fits')
    #
    #first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    #first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    ## 2nd order polynomial
    # file_out=file_in.replace('.fits','_base2.fits')
    #file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    #first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    #print("Now NH3(2,2)")
    #a_rms = [  0, 260, 520, 730]
    #b_rms = [150, 380, 610, 850]
    #index_rms=first_look.create_index( a_rms, b_rms)
    #index_peak=np.arange(380,520)
    #line='NH3_22'
    #file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    #file_out=file_in.replace(file_extension+'.fits',
    #                             '_base'+file_extension+'.fits')
    #first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    #first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    ## 2nd order polynomial
    #file_out=file_in.replace('.fits','_base2.fits')
    #file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    #first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    #print("Now NH3(3,3)")
    #a_rms = [ 10, 250, 530]
    #b_rms = [210, 310, 930]
    #index_rms=first_look.create_index( a_rms, b_rms)
    #index_peak=np.arange(410,540)
    #line='NH3_33'
    #file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    #file_out=file_in.replace(file_extension+'.fits',
    #                             '_base'+file_extension+'.fits')
    #first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    #first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    ##2nd order polynomial
    #file_out=file_in.replace('.fits','_base2.fits')
    #file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    #first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    #print("Now CCS")
    #a_rms = [  0, 260]
    #b_rms = [200, 490]
    #index_rms=first_look.create_index( a_rms, b_rms)
    #index_peak=np.arange(220,250)
    #line='C2S'
    #file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    #file_out=file_in.replace(file_extension+'.fits',
    #                             '_base'+file_extension+'.fits')
    #first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    #first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    #print("Now HC5N")
    # HC5N channel range must be updated
    #a_rms = [  0, 500]
    #b_rms = [380, 545]
    #index_rms=first_look.create_index( a_rms, b_rms)
    #index_peak=np.arange(400,480)
    #line='HC5N'
    #file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    #file_out=file_in.replace(file_extension+'.fits',
    #                             '_base'+file_extension+'.fits')
    #first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    #first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    #print("Now HC7N 21-20")
    # HC7N channel range must be updated
    #a_rms = [  0, 160, 480]
    #b_rms = [115, 360, 525]
    #index_rms=first_look.create_index( a_rms, b_rms)
    #index_peak=np.arange(400,460)
    #line='HC7N_21_20'
    #file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    #file_out=file_in.replace(file_extension+'.fits',
    #                             '_base'+file_extension+'.fits')
    #first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    #first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    #print("Now HC7N 22-21")
    # HC7N channel range must be updated
    #a_rms = [  0, 480]
    #b_rms = [360, 525]
    #index_rms=first_look.create_index( a_rms, b_rms)
    #index_peak=np.arange(400,460)
    #line='HC7N_22_21'
    #file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    #file_out=file_in.replace(file_extension+'.fits',
    #                             '_base'+file_extension+'.fits')
    #first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    #first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)

def FirstLook_B18(file_extension='_all'):
    """
    Function to create First Look products for B18. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name='B18'
    vsys = 6.*u.km/u.s
    throw = 2.0*u.km/u.s

    print("Now NH3(1,1)")
    
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)

    s = SpectralCube.read(file_in)
    s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
    spaxis = s.spectral_axis.value
    index_rms = baseline.ammoniaWindow(spaxis,spaxis,window=4,v0=vsys.value)
    index_peak= ~baseline.tightWindow(spaxis,spaxis,window=3,v0=vsys.value)
    first_look.peak_rms( file_in, index_rms=index_rms, index_peak=index_peak)

    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']

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
        first_look.peak_rms( file_in, index_rms=index_rms, 
                             index_peak=index_peak)

    #region_name='B18'
    #print("Now NH3(1,1)")
    #a_rms = [  0, 115, 280, 385, 490, 655]
    #b_rms = [ 80, 230, 345, 455, 625, 760]
    #index_rms=first_look.create_index( a_rms, b_rms)
    #index_peak=np.arange(352,381)
    #file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    #file_out=file_in.replace(file_extension+'.fits',
    #                         '_base'+file_extension+'.fits')
    #first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    #first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    #print("Now NH3(2,2)")
    #a_rms = [   0, 440]
    #b_rms = [ 409, 870]
    #index_rms=first_look.create_index( a_rms, b_rms)
    #index_peak=np.arange(420,435)
    #line='NH3_22'
    #file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    #file_out=file_in.replace(file_extension+'.fits',
    #                             '_base'+file_extension+'.fits')
    #first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    #first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    #print("Now NH3(3,3)")
    #a_rms = [   0, 530]
    #b_rms = [ 409, 960]
    #index_rms=first_look.create_index( a_rms, b_rms)
    #index_peak=np.arange(410,485)
    #line='NH3_33'
    #file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    #file_out=file_in.replace(file_extension+'.fits',
    #                             '_base'+file_extension+'.fits')
    #first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    #first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    #print("Now CCS")
    #a_rms = [   0, 245]
    #b_rms = [ 210, 490]
    #index_rms=first_look.create_index( a_rms, b_rms)
    #index_peak=np.arange(225,243)
    #line='C2S'
    #file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    #file_out=file_in.replace(file_extension+'.fits',
    #                             '_base'+file_extension+'.fits')
    #first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    #first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    #print("Now HC5N")
    #a_rms = [  10, 435]
    #b_rms = [ 409, 540]
    #index_rms=first_look.create_index( a_rms, b_rms)
    #index_peak=np.arange(414,430)
    #line='HC5N'
    #file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    #file_out=file_in.replace(file_extension+'.fits',
    #                             '_base'+file_extension+'.fits')
    #first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    #first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    #print("Now HC7N_21_20")
    #a_rms = [  10, 435]
    #b_rms = [ 409, 540]
    #index_rms=first_look.create_index( a_rms, b_rms)
    #index_peak=np.arange(412,430)
    #line='HC7N_21_20'
    #file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    #file_out=file_in.replace(file_extension+'.fits',
    #                             '_base'+file_extension+'.fits')
    #first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    #first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    #
    #print("Now HC7N_22_21")
    #a_rms = [  10, 435]
    #b_rms = [ 409, 540]
    #index_rms=first_look.create_index( a_rms, b_rms)
    #index_peak=np.arange(412,430)
    #line='HC7N_22_21'
    #file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    #file_out=file_in.replace(file_extension+'.fits',
    #                             '_base'+file_extension+'.fits')
    #first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    #first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)


def FirstLook_L1688(file_extension='_all'):
    """
    Function to create First Look products for L1688. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name='L1688'
    vsys = 3.5*u.km/u.s
    throw = 5*u.km/u.s

    print("Now NH3(1,1)")
    #a_rms = [  0, 121, 290, 404, 505, 665]
    #b_rms = [ 74, 239, 332, 447, 611, 749]
    #index_rms=first_look.create_index( a_rms, b_rms)
    #index_peak=np.arange(350,377)
    
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)
    # file_out=file_in.replace(file_extension+'.fits',
    #                          '_base'+file_extension+'.fits')
    # first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    s = SpectralCube.read(file_in)
    s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
    spaxis = s.spectral_axis.value
    index_rms = baseline.ammoniaWindow(spaxis,spaxis,window=4,v0=vsys.value)
    index_peak= ~baseline.tightWindow(spaxis,spaxis,window=3,v0=vsys.value)
    first_look.peak_rms( file_in, index_rms=index_rms, index_peak=index_peak)

    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']


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

        #file_out=file_in.replace(file_extension+'.fits',
        #                         '_base'+file_extension+'.fits')
        # first_look.baseline( file_in, file_out, 
        #                               index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_in, index_rms=index_rms, 
                             index_peak=index_peak)


    # print("Now NH3(2,2)")
    # a_rms = [   0, 349]
    # b_rms = [ 285, 649]
    # index_rms=first_look.create_index( a_rms, b_rms)
    # index_peak=np.arange(298,342)
    # line='NH3_22'
    # file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    # file_out=file_in.replace(file_extension+'.fits',
    #                              '_base'+file_extension+'.fits')
    # first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    # first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    # #
    # print("Now NH3(3,3)")
    # a_rms = [   0, 395]
    # b_rms = [ 272, 649]
    # index_rms=first_look.create_index( a_rms, b_rms)
    # index_peak=np.arange(298,342)
    # line='NH3_33'
    # file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    # file_out=file_in.replace(file_extension+'.fits',
    #                              '_base'+file_extension+'.fits')
    # first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    # first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    # #
    # print("Now CCS")
    # a_rms = [   0, 369]
    # b_rms = [ 278, 649]
    # index_rms=first_look.create_index( a_rms, b_rms)
    # index_peak=np.arange(307,325)
    # line='C2S'
    # file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    # file_out=file_in.replace(file_extension+'.fits',
    #                              '_base'+file_extension+'.fits')
    # first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    # first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    # #
    # print("Now HC5N")
    # a_rms = [   0, 358]
    # b_rms = [ 288, 649]
    # index_rms=first_look.create_index( a_rms, b_rms)
    # index_peak=np.arange(306,317)
    # line='HC5N'
    # file_in = '{0}/{0}_{1}{2}.fits'.format(region_name,line,file_extension)
    # file_out=file_in.replace(file_extension+'.fits',
    #                              '_base'+file_extension+'.fits')
    # first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    # first_look.peak_rms( file_out, index_rms=index_rms, index_peak=index_peak)
    # #
    # #HC7N (21-20) shows an absorption feature at ~ 91 km/s (at 23.6951 GHz)
    # #from its rest frequency (used 23.6879 GHz). There's no emission line.
    # #Below are the channel indeces for the absorption feature.
    # #a_rms = [  0, 520]
    # #b_rms = [480, 650]
    # #index_peak = np.arange(485,510)
    # #
    # #The code didn't produce the fits file for HC7N (22-21).

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
    vsys = 7.9*u.km/u.s
    throw = 2.0*u.km/u.s

    print("Now NH3(1,1)")
    
    file_in='{0}/{0}_NH3_11{1}.fits'.format(region_name,file_extension)

    s = SpectralCube.read(file_in)
    s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
    spaxis = s.spectral_axis.value
    index_rms = baseline.ammoniaWindow(spaxis,spaxis,window=4,v0=vsys.value)
    index_peak= ~baseline.tightWindow(spaxis,spaxis,window=3,v0=vsys.value)
    first_look.peak_rms( file_in, index_rms=index_rms, index_peak=index_peak)

    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']

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
        first_look.peak_rms( file_in, index_rms=index_rms, 
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
        
def FirstLook_Cepheus_L1228(file_extension='_all'):
    """
    Function to create First Look products for Cepheus L1228. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name = 'Cepheus_L1228'
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
    vsys = -8.0*u.km/u.s
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

def FirstLook_Cepheus_L1251(file_extension='_all'):
    """
    Function to create First Look products for Cepheus_L1251. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name = 'Cepheus_L1251'
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

def FirstLook_IC5146(file_extension='_all'):
    """
    Function to create First Look products for IC5146. The file_extension 
    parameter is used to select the proper files to be processed. 
    """
    region_name = 'IC5146'
    print("Now NH3(1,1)")
    a_rms = [  0, 135, 290, 405, 505, 660]
    b_rms = [ 70, 235, 340, 445, 615, 740]
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
    vsys = 4.0*u.km/u.s
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

def FirstLook_SerMWC(file_extension='_all'):
    """
    Function to create First Look products for Serpens_Aquila. The 
    file_extension parameter is used to select the proper files to be 
    processed. 
    """
    region_name='Serpens_MWC297'
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

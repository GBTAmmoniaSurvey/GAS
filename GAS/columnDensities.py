import astropy.io.fits as fits
import numpy as np
import os
import astropy.constants as c
import astropy.units as u
import glob
from scipy.ndimage import binary_opening
import aplpy

import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

from config import plottingDictionary

#####################################################
# Calculate column densities of non-NH3 lines in GAS
# based on results of Gaussian fit analysis
# To do: possibly include hyperfine analysis of HC5N
#####################################################

def calc_jt(temp, nu):
    factor = np.exp(c.h * nu/(c.k_B * temp))-1.
    jt = c.h * nu / c.k_B  / factor
    return jt

def calc_tau(tmb,tex, nu):
    tbg = 2.73*u.K
    factor = tmb/(calc_jt(tex,nu)-calc_jt(tbg,nu))
    tau = (-1.)*np.log(1. - factor)
    return tau

def calc_q_linear(tex,b):
    # Calc Q
    # Qs seem reasonable for HC5N compared with Splatalogue values, poorer agreement with CCS
    j_arr = np.array(range(100))
    q=0
    for j in j_arr:
        q += (2.*j + 1)*np.exp(-1.*c.h.cgs*b.cgs*j*(j+1)/(c.k_B.cgs*tex))
    return q

def calc_q_ccs(tex):
    # For each J have three different states (N = J +/- 1)
    ccs_file = os.path.join(os.path.dirname(__file__),'CCS_lvls.txt')
    Eu, Nu, Ju = np.loadtxt(ccs_file,unpack=True)
    gu = 2.*Ju + 1.
    q = np.sum(gu * np.exp((-1.)*Eu/tex.value))
    return q

def calc_n_mangum(tmb,tex,dv,params,q):
    nu = params['nu']
    tau = calc_tau(tmb,tex,nu)
    jup = params['jup']
    Eup = params['Eup']
    mu = params['mu']
    Sij = params['Sij']
    b = params['b']
    tbg = 2.73*u.K
    #q = calc_q_linear(tex,b)
    n = tau * 3.*c.h/(8.*np.pi**3.*Sij*mu**2)*q/(2.*jup+1)*\
        np.exp(Eup/tex)*(np.exp(c.h*nu/(c.k_B*tex))-1.)**(-1.)*\
        dv*np.sqrt(np.pi/(4.*np.log(2.)))
    return n, tau

def calc_n_mom0(mom0,tmb,tex,nu,tex_arr,n_arr):
    # Use assumed excitation temperature and Yancy's text files to convert from integrated intensity
    # to column density (read tex_arr, n_arr based on line)
    tau = calc_tau(tmb,tex,nu)
    itex = np.where(tex_arr == tex.value) # Can do this if tex is less than 2 decimal precision
    n = mom0 * n_arr[itex][0]
    # Tau correction where needed
    corr = np.where(tau > 0.1)
    n[corr] = n[corr] * tau[corr]/(1.-np.exp(-(1.)*tau[corr]))
    return n

def calc_column_densities_fits(region='B18',file_extension='all_rebase3',tex=7.*u.K):
    if file_extension:
        root = file_extension
    else:
        root = 'all'
    
    lines = ['HC5N','C2S','HC7N_21_20']
    hc5n_params = {'nu' : 23.9638968*u.GHz,
                   'jlow' : 8,
                   'jup' : 9,
                   'Elow' : 4.6003*u.K,
                   'Eup' : 5.7505*u.K,
                   'mu' : 4.33e-18*u.esu*u.cm,
                   'Sij' : 0.47, 
                   'b' : 1331.33*u.MHz}
    c2s_params = {'nu' : 22.3440308*u.GHz,
                  'jlow' : 1, 
                  'jup' : 2,
                  'Elow' : 0.5336*u.K,
                  'Eup' : 1.857*u.K,
                  'mu' : 2.88e-18*u.esu*u.cm,
                  'Sij' : 0.4,
                  'b' : 6477.75*u.MHz}
    hc7n_21_params = {'nu': 23.6878974*u.GHz,
                  'jlow' : 20, 
                  'jup' : 21,
                  'Elow' : 11.3684*u.K,
                  'Eup' : 12.5052*u.K,
                  'mu' : 5.0e-18*u.esu*u.cm,
                  'Sij' : 0.489, 
                  'b' : 564.0*u.MHz}
    line_params = {'HC5N' : hc5n_params, 'C2S':c2s_params, 'HC7N_21_20':hc7n_21_params}

    for line in lines:
        gparamfits = '{0}/{0}_{2}_{1}_param_cube_masked.fits'.format(region,root,line)
        mom0file = '{0}/{0}_{1}_{2}_mom0_QA.fits'.format(region,line,root)
        colfile = '{0}/parameterMaps/{0}_{1}_{2}_N_masked.fits'.format(region,line,root)
        taufile = '{0}/parameterMaps/{0}_{1}_{2}_tau_masked.fits'.format(region,line,root)
        # Make sure files exist
        if os.path.isfile(gparamfits):
            gparam_hdu = fits.open(gparamfits)
            gparam_data = gparam_hdu[0].data
            mom0_hdu = fits.open(mom0file)
            header = mom0_hdu[0].header
            mom0 = mom0_hdu[0].data
            tmb_fit = gparam_data[0] * u.K
            sigv_fit = gparam_data[2] * u.km/u.s
            fwhm = 2.*np.sqrt(2.*np.log(2.))*sigv_fit
            params = line_params[line]
            if line == 'C2S':
                q = calc_q_ccs(tex)
            else:
                q = calc_q_linear(tex,params['b'])
            ncol, tau = calc_n_mangum(tmb_fit,tex,fwhm,params,q)
            # Edit header
            #rm_key=['CDELT3', 'CUNIT3', 'CTYPE3', 'CRVAL3']
            #for key_i in rm_key:
            #    header.remove(key_i)
            header['NAXIS'] = 2
            header['WCSAXES'] = 2
            header['BUNIT'] = 'cm-2'
            # Write out files
            new_hdu = fits.PrimaryHDU(ncol.cgs.value,header=header)
            new_hdu.writeto(colfile,overwrite=True)        
            header['BUNIT'] = ''
            new_hdu2 = fits.PrimaryHDU(tau.value,header=header)
            new_hdu2.writeto(taufile,overwrite=True)


def calc_all_columns(file_extension='all_rebase3',tex=7.*u.K,release='all'):
    RegionCatalog = catalogs.GenerateRegions(release=release)

    for ThisRegion in RegionCatalog:
        region = ThisRegion['Region name']
        calc_column_densities_fits(region=region,file_extension=file_extension,tex=tex)


def plot_property_maps(regions=None,file_extension='all_rebase3',release='all'):
    # Get list of regions - run from images/ directory
    # Assume directories correspond to regions to be imaged
    # Update - use region list
    if regions is None:
        RegionCatalog = catalogs.GenerateRegions(release=release)
    else:
        RegionCatalog = catalogs.GenerateRegions(release=release)
        keep = [idx for idx, row in enumerate(RegionCatalog) if row['Region name'] in regions]
        RegionCatalog = RegionCatalog[keep]

    ext_list  = [0,1,2]
    label_list = ['$T_B (K)','$v_\mathrm{LSR}$ (km s$^{-1}$)','$\sigma_v$ (km s$^{-1}$)']
    file_list = ['T_B','vlsr','sigv']
    ctable_list = ['plasma','RdYlBu_r','plasma'] #'YlGnBu_r'
    text_color='black'
    text_size = 12
    beam_color='#d95f02'  # previously used '#E31A1C'
    # Try single set of contours for first look images
    w11_step = 0.4
    cont_levs=2**np.arange( 0,20)*w11_step
    w11_lw   = 0.5
    # Masking of small (noisy) regions
    selem = np.array([[0,1,0],[1,1,1],[0,1,0]])

    line_list = ['HC5N','C2S','HC7N_21_20','HC7N_22_21','NH3_33']

    for ThisRegion in RegionCatalog:
        region = ThisRegion['Region name']
        if os.path.isdir(region):
            print region
            # Use NH3 (1,1) moment maps for contours? 
            file_w11='{0}/{0}_NH3_11_{1}_mom0_QA_trim.fits'.format(region,file_extension)
            plot_param = plottingDictionary[region]
            for line in line_list:
                gparamfits = '{0}/{0}_{2}_{1}_param_cube_masked.fits'.format(region,file_extension,line)
                if os.path.isfile(gparamfits):
                    par_hdu = fits.open(gparamfits)
                    # Get NH3 (1,1) moment contours
                    LowestContour= cont_levs[0]*0.5
                    nh3mom0_hdu = fits.open(file_w11)
                    nh3mom0 = nh3mom0_hdu[0].data
                    mask = binary_opening(nh3mom0 > LowestContour, selem)
                    MaskedMap = mask*nh3mom0
                    nh3mom0_hdu[0].data = MaskedMap
                    for i in range(len(ext_list)):
                        # Use percentiles to set initial plot colourscale ranges
                        # Need to omit zeros from calculation
                        plane = par_hdu[0].data[i,:,:]
                        v_min=np.nanpercentile(plane[np.where(plane !=0)],2.5)
                        v_max=np.nanpercentile(plane[np.where(plane !=0)],97.5)
                        fig=aplpy.FITSFigure(par_hdu,slices=[i])
                        fig.show_colorscale(cmap=ctable_list[i],vmin=v_min, vmax=v_max)
                        # For some reason having set_nan_color *after* colorbar messes up the tick locations! 
                        fig.set_nan_color('0.99')
                        # add colorbar
                        fig.add_colorbar()
                        fig.colorbar.show(box_orientation='horizontal', width=0.1, pad=0.0,
                                          location='top')
                        fig.colorbar.set_font(family='sans_serif',size=text_size)
                        fig.colorbar.set_axis_label_font(family='sans_serif',size=text_size)
                        #
                        fig.show_contour(nh3mom0_hdu, colors='gray', levels=cont_levs, linewidths=w11_lw)
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
                        # Scale bar
                        # magic line of code to obtain scale in arcsec obtained from
                        # http://www.astropy.org/astropy-tutorials/Quantities.html
                        ang_sep =  (plot_param['scalebar_size'].to(u.au)/plot_param['distance']).to(u.arcsec, equivalencies=u.dimensionless_angles())
                        fig.add_scalebar(ang_sep.to(u.degree))
                        fig.scalebar.set_corner(plot_param['scalebar_pos'])
                        fig.scalebar.set_font(family='sans_serif',size=text_size)
                        fig.scalebar.set(color=text_color)
                        fig.scalebar.set_label('{0:4.2f}'.format(plot_param['scalebar_size']))
                        # Labels
                        label = label_list[i]
                        label_loc = plot_param['label_loc']
                        label_ha = plot_param['label_ha']
                        fig.add_label(label_loc[0],label_loc[1],
                                      '{0}\n{1}'.format(region,label),
                                      relative=True, color=text_color,
                                      horizontalalignment=label_ha,
                                      family='sans_serif',size=text_size)
                        # fig.set_system_latex(True)
                        fig.save('figures/{0}_{1}_{2}_{3}.pdf'.format(region,line,file_extension,file_list[i]),
                                 adjust_bbox=True,dpi=200)
                        fig.close()
                else:
                    print('File {0} not found'.format(gparamfits))

    # Next plot column density maps
    for region in region_list:
        if os.path.isdir(region):
            print region
            # Use NH3 (1,1) moment maps for contours? 
            file_w11='{0}/{0}_NH3_11_{1}_mom0_QA_trim.fits'.format(region,file_extension)
            plot_param = plottingDictionary[region]
            for line in line_list:
                colfits = '{0}/parameterMaps/{0}_{2}_{1}_N_masked.fits'.format(region,file_extension,line)
                if os.path.isfile(colfits):
                    col_hdu = fits.open(colfits)
                    # Get NH3 (1,1) moment contours
                    LowestContour= cont_levs[0]*0.5
                    nh3mom0_hdu = fits.open(file_w11)
                    nh3mom0 = nh3mom0_hdu[0].data
                    mask = binary_opening(nh3mom0 > LowestContour, selem)
                    MaskedMap = mask*nh3mom0
                    nh3mom0_hdu[0].data = MaskedMap
                    # Plot log of column for clarity
                    log_data = np.log10(col_hdu[0].data)
                    col_hdu[0].data = log_data
                    # Use percentiles to set initial plot colourscale ranges
                    v_min=np.nanpercentile(col_hdu[0].data[np.where(col_hdu[0].data > 0)],0.5)
                    v_max=np.nanpercentile(col_hdu[0].data[np.where(col_hdu[0].data > 0)],99.5)
                    v_mid = 0
                    fig=aplpy.FITSFigure(col_hdu)
                    fig.show_colorscale(cmap='Blues',vmin=v_min,vmax=v_max)
                    # For some reason having set_nan_color *after* colorbar messes up the tick locations! 
                    fig.set_nan_color('0.95')
                    # add colorbar
                    fig.add_colorbar()
                    fig.colorbar.show(box_orientation='horizontal', width=0.1, pad=0.0,
                                      location='top')
                    fig.colorbar.set_font(family='sans_serif',size=text_size)
                    fig.colorbar.set_axis_label_font(family='sans_serif',size=text_size)
                    #
                    fig.show_contour(nh3mom0_hdu, colors='gray', levels=cont_levs, linewidths=w11_lw)
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
                    # Scale bar
                    # magic line of code to obtain scale in arcsec obtained from
                    # http://www.astropy.org/astropy-tutorials/Quantities.html
                    ang_sep =  (plot_param['scalebar_size'].to(u.au)/plot_param['distance']).to(u.arcsec, equivalencies=u.dimensionless_angles())
                    fig.add_scalebar(ang_sep.to(u.degree))
                    fig.scalebar.set_corner(plot_param['scalebar_pos'])
                    fig.scalebar.set_font(family='sans_serif',size=text_size)
                    fig.scalebar.set(color=text_color)
                    fig.scalebar.set_label('{0:4.2f}'.format(plot_param['scalebar_size']))
                    # Labels
                    label_loc = plot_param['label_loc']
                    label_ha = plot_param['label_ha']
                    fig.add_label(label_loc[0],label_loc[1],
                                  '{0}\nlog N({1})'.format(region,line),
                                  relative=True, color=text_color,
                                  horizontalalignment=label_ha,
                                  family='sans_serif',size=text_size)
                    # fig.set_system_latex(True)
                    fig.save('figures/{0}_{1}_{2}_N.pdf'.format(region,line,file_extension),
                             adjust_bbox=True,dpi=200)
                    fig.close()
                else:
                    print('File {0} not found'.format(gparamfits))

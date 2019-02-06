import astropy.io.fits as fits
import numpy as np
import os
import astropy.constants as c
import astropy.units as u
import glob

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

def calc_column_densities_fits(region='B18',file_extension=None,tex=7.*u.K):
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
        colfile = '{0}/{0}_{1}_{2}_N.fits'.format(region,line,root)
        taufile = '{0}/{0}_{1}_{2}_tau.fits'.format(region,line,root)
        # Make sure files exist
        if os.path.isfile(gparamfits):
            gparam_hdu = fits.open(gparamfits)
            header = gparam_hdu[0].header
            gparam_data = gparam_hdu[0].data
            mom0_hdu = fits.open(mom0file)
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
            header['BUNIT'] = 'cm-2'
            # Write out files
            new_hdu = fits.PrimaryHDU(ncol.cgs.value,header=header)
            new_hdu.writeto(colfile,overwrite=True)        
            header['BUNIT'] = ''
            new_hdu2 = fits.PrimaryHDU(tau.value,header=header)
            new_hdu2.writeto(taufile,overwrite=True)


def calc_all_columns(file_extension='all_rebase3',tex=7.*u.K):
    # Get region list from directories in GAS/imaging/
    region_list = glob.glob("*/")
    for i in range(len(region_list)):
        region_list[i] = region_list[i].strip("/")
    if 'figures' in region_list: region_list.remove('figures')

    for region in region_list:
        calc_column_densities_fits(region=region,file_extension=file_extension,tex=tex)

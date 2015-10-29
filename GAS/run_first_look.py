from . import first_look
import numpy as np
from spectral_cube import SpectralCube
import astropy.units as u
def FirstLook_OrionA():
    print("Now NH3(1,1)")
    a_rms = [  0, 158, 315, 428, 530, 693]
    b_rms = [ 60, 230, 327, 438, 604, 735]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(326,430)
    file_in='OrionA/OrionA_NH3_11.fits'
    # 1st order polynomial
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
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
    file_in='OrionA/OrionA_NH3_22.fits'
    # 1st order polynomial
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
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
    file_in='OrionA/OrionA_NH3_33.fits'
    # 1st order polynomial
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
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
    file_in='OrionA/OrionA_C2S.fits'
    file_out=file_in.replace('.fits','_base2.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC5N")
    # HC5N channel range must be updated
    a_rms = [  0, 500]
    b_rms = [380, 545]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(400,480)
    file_in='OrionA/OrionA_HC5N.fits'
    file_out=file_in.replace('.fits','_base2.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC7N 21-20")
    # HC7N channel range must be updated
    a_rms = [  0, 160, 480]
    b_rms = [115, 360, 525]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(400,460)
    file_in='OrionA/OrionA_HC7N_21_20.fits'
    file_out=file_in.replace('.fits','_base2.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC7N 22-21")
    # HC7N channel range must be updated
    a_rms = [  0, 480]
    b_rms = [360, 525]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(400,460)
    file_in='OrionA/OrionA_HC7N_22_21.fits'
    file_out=file_in.replace('.fits','_base2.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

def FirstLook_B18():
    print("Now NH3(1,1)")
    a_rms = [  0, 115, 280, 385, 490, 655]
    b_rms = [ 80, 230, 345, 455, 625, 760]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(352,381)
    file_in='B18/B18_NH3_11.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #

    print("Now NH3(2,2)")
    a_rms = [   0, 440]
    b_rms = [ 409, 870]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(420,435)
    file_in='B18/B18_NH3_22.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now NH3(3,3)")
    a_rms = [   0, 530]
    b_rms = [ 409, 960]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(410,485)
    file_in='B18/B18_NH3_33.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now CCS")
    a_rms = [   0, 245]
    b_rms = [ 210, 490]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(225,243)
    file_in='B18/B18_C2S.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC5N")
    a_rms = [  10, 435]
    b_rms = [ 409, 540]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(414,430)
    file_in='B18/B18_HC5N.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC7N_21_20")
    a_rms = [  10, 435]
    b_rms = [ 409, 540]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(412,430)
    file_in='B18/B18_HC7N_21_20.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC7N_22_21")
    a_rms = [  10, 435]
    b_rms = [ 409, 540]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(412,430)
    file_in='B18/B18_HC7N_22_21.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)


def FirstLook_L1688():
    print("Now NH3(1,1)")
    a_rms = [  0, 121, 290, 404, 505, 665]
    b_rms = [ 74, 239, 332, 447, 611, 749]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(350,377)
    file_in='L1688/L1688_NH3_11.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now NH3(2,2)")
    a_rms = [   0, 349]
    b_rms = [ 285, 649]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(298,342)
    file_in='L1688/L1688_NH3_22.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now NH3(3,3)")
    a_rms = [   0, 395]
    b_rms = [ 272, 649]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(298,342)
    file_in='L1688/L1688_NH3_33.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now CCS")
    a_rms = [   0, 369]
    b_rms = [ 278, 649]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(307,325)
    file_in='L1688/L1688_C2S.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC5N")
    a_rms = [   0, 358]
    b_rms = [ 288, 649]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(306,317)
    file_in='L1688/L1688_HC5N.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    #HC7N (21-20) shows an absorption feature at ~ 91 km/s (at 23.6951 GHz)
    #from its rest frequency (used 23.6879 GHz). There's no emission line.
    #Below are the channel indeces for the absorption feature.
    #a_rms = [  0, 520]
    #b_rms = [480, 650]
    #index_peak = np.arange(485,510)
    #
    #The code didn't produce the fits file for HC7N (22-21).

def FirstLook_SerAqu():
    print("Now NH3(1,1)")
    a_rms = [  0, 150, 310, 420, 530, 690]
    b_rms = [ 60, 230, 330, 440, 610, 780]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(340,420)
    file_in='Serpens_Aquila/Serpens_Aquila_NH3_11.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now NH3(2,2)")
    a_rms = [  0, 230, 460, 665]
    b_rms = [150, 380, 610, 820]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(380,460)
    file_in='Serpens_Aquila/Serpens_Aquila_NH3_22.fits'
    # 1st order polynomial
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now NH3(3,3)")
    a_rms = [  0, 300]
    b_rms = [220, 490]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(380,460)
    file_in='Serpens_Aquila/Serpens_Aquila_NH3_33.fits'
    # 1st order polynomial
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now CCS")
    a_rms = [  0, 260]
    b_rms = [200, 490]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(220,250)
    file_in='Serpens_Aquila/Serpens_Aquila_C2S.fits'
    # 1st order polynomial
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC5N")
    a_rms = [  0, 260]
    b_rms = [200, 490]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(220,250)
    file_in='Serpens_Aquila/Serpens_Aquila_HC5N.fits'
    # 1st order polynomial
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

def FirstLook_L1455():
    print("Now NH3(1,1)")
    a_rms = [   0, 140, 300, 410, 520, 680]
    b_rms = [ 105, 270, 370, 480, 630, 745]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(350,430)
    file_in='L1455/L1455_NH3_11.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

    print("Now NH3(2,2)")
    a_rms = [   0, 340]
    b_rms = [ 290, 648]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(260,400)
    file_in='L1455/L1455_NH3_22.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

    print("Now NH3(3,3)")
    a_rms = [   0, 340]  # No lines. Using the same as NH3(2,2)
    b_rms = [ 290, 648]  # No lines. Using the same as NH3(2,2)
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(260,400)  # No lines. Using the same as NH3(2,2)
    file_in='L1455/L1455_NH3_33.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

    print("Now CCS")
    a_rms = [   0, 350]  
    b_rms = [ 290, 648]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(309,334)  
    file_in='L1455/L1455_C2S.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

    print("Now HC5N")
    a_rms = [   0, 350]  
    b_rms = [ 290, 648]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(315,325)  
    file_in='L1455/L1455_HC5N.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

    print("Now HC7N_21_20")
    a_rms = [   0, 180]  
    b_rms = [ 130, 275]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(128,147)  
    file_in='L1455/L1455_HC7N_21_20.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

    print("Now HC7N_22_21")
    a_rms = [   0, 340]  # No lines. Using the same as HC7N_21_20
    b_rms = [ 290, 648]  # No lines. Using the same as HC7N_21_20
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(308,328)  # No lines. Using the same as HC7N_21_20
    file_in='L1455/L1455_HC7N_22_21.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)


def FirstLook_NGC1333():
    print("Now NH3(1,1)")
    a_rms = [  0, 158, 315, 428, 530, 693, 751]
    b_rms = [ 60, 230, 327, 438, 604, 735, 760]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(326,430)
    file_in='NGC1333/NGC1333_NH3_11.fits'
    # 1st order polynomial
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    print("Now NH3(2,2)")
    a_rms = [  0, 190, 360, 600]
    b_rms = [70, 300, 470, 640]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(380,520)
    # file_in='NGC1333/NGC1333_NH3_22.fits'
    # # 1st order polynomial
    # file_out=file_in.replace('.fits','_base1.fits')
    # file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    # first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    ## 2nd order polynomial
    #file_out=file_in.replace('.fits','_base2.fits')
    #file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    #first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    # print("Now NH3(3,3)")
    # a_rms = [ 10, 190, 420]
    # b_rms = [70, 360, 500]
    # index_rms=first_look.create_index( a_rms, b_rms)
    # index_peak=np.arange(410,540)
    # file_in='NGC1333/NGC1333_NH3_33.fits'
    # 1st order polynomial
    # file_out=file_in.replace('.fits','_base1.fits')
    # file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    # first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 8.5*u.km/u.s
    throw = 8*u.km/u.s
    for line in linelist:
        file_in = 'NGC1333/NGC1333_{0}.fits'.format(line)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+2*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-2*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace('.fits','_base1.fits')
        file_new=first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_new, index_rms=index_rms, 
                             index_peak=index_peak)
        

def FirstLook_B1():
    print("Now NH3(1,1)")
    a_rms = [  0, 130, 290, 400, 500, 660]
    b_rms = [ 70, 240, 340, 440, 620, 740]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(340,400)
    file_in='B1/B1_NH3_11.fits'
    # 1st order polynomial
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    print("Now NH3(2,2)")
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 6.6*u.km/u.s
    throw = 2.0*u.km/u.s
    for line in linelist:
        file_in = 'B1/B1_{0}.fits'.format(line)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+3*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-3*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace('.fits','_base1.fits')
        file_new=first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_new, index_rms=index_rms, 
                             index_peak=index_peak)
        
def FirstLook_IC348():
    print("Now NH3(1,1)")
    a_rms = [  0, 130, 290, 400, 500, 660]
    b_rms = [ 70, 240, 340, 440, 620, 740]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(340,400)
    file_in='IC348/IC348_NH3_11.fits'
    # 1st order polynomial
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    print("Now NH3(2,2)")
    linelist = ['NH3_22','NH3_33','C2S','HC5N','HC7N_21_20','HC7N_22_21']
    vsys = 9.0*u.km/u.s
    throw = 2.0*u.km/u.s
    for line in linelist:
        file_in = 'IC348/IC348_{0}.fits'.format(line)
        s = SpectralCube.read(file_in)
        s = s.with_spectral_unit(u.km/u.s,velocity_convention='radio')
        a_rms = [s.closest_spectral_channel(vsys+3*throw),
                 s.closest_spectral_channel(vsys-throw)]
        b_rms = [s.closest_spectral_channel(vsys+throw),
                 s.closest_spectral_channel(vsys-3*throw)]
        index_peak = np.arange(s.closest_spectral_channel(vsys+3*u.km/u.s),
                              s.closest_spectral_channel(vsys-3*u.km/u.s))
        index_rms=first_look.create_index( a_rms, b_rms)

        file_out=file_in.replace('.fits','_base1.fits')
        file_new=first_look.baseline( file_in, file_out, 
                                      index_clean=index_rms, polyorder=1)
        first_look.peak_rms( file_new, index_rms=index_rms, 
                             index_peak=index_peak)
        

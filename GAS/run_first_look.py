from . import first_look
import numpy as np


def FirstLook_OrionA():
    print("Now NH3(1,1)")
    a_rms = [ 50, 300, 440, 555, 670, 820, 875]
    b_rms = [170, 350, 450, 565, 720, 860, 910]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(445,560)
    file_in='OrionA/OrionA_NH3_11.fits'
    # 1st order polynomial
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    # 2nd order polynomial
    file_out=file_in.replace('.fits','_base2.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    print("Now NH3(2,2)")
    a_rms = [  0, 270, 575, 745]
    b_rms = [225, 420, 660, 935]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(425,562)
    file_in='OrionA/OrionA_NH3_22.fits'
    # 1st order polynomial
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    # 2nd order polynomial
    file_out=file_in.replace('.fits','_base2.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    print("Now NH3(3,3)")
    a_rms = [  0, 290, 425, 600]
    b_rms = [255, 410, 445, 945]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(440,580)
    file_in='OrionA/OrionA_NH3_33.fits'
    # 1st order polynomial
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #2nd order polynomial
    file_out=file_in.replace('.fits','_base2.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    print("Now HC5N")
    a_rms = [  0, 550]
    b_rms = [420, 945]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(450,560)
    file_in='OrionA/OrionA_HC5N.fits'
    # 1st order polynomial
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    # 2nd order polynomial
    file_out=file_in.replace('.fits','_base2.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=2)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)


def FirstLook_NGC1333():
    print("Now NH3(1,1)")
    a_rms = [  0, 165, 239, 325, 435, 540, 700]
    b_rms = [100, 200, 270, 370, 475, 630, 795]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(370,450)
    file_in='NGC1333/NGC1333_NH3_11.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

def FirstLook_B18():
    print("Now NH3(1,1)")
    a_rms = [  0, 145, 230, 310, 420, 520, 690]
    b_rms = [100, 205, 270, 375, 475, 640, 795]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(370,415)
    file_in='B18/B18_NH3_11.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)


def FirstLook_L1688():
    print("Now NH3(1,1)")
    a_rms = [  0, 145, 230, 310, 420, 525, 690]
    b_rms = [ 95, 210, 265, 360, 470, 640, 795]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(370,415)
    file_in='L1688/L1688_NH3_11.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

def FirstLook_SerAqu():
    print("Now NH3(1,1)")
    a_rms = [  0, 150, 310, 420, 530, 690]
    b_rms = [ 60, 230, 330, 440, 610, 780]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(340,420)
    file_in='Serpens_Aquila_test/Serpens_Aquila_NH3_11.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

def FirstLook_L1544():
    print("Now NH3(1,1)")
    a_rms = [  0,  55, 135, 185, 235, 310]
    b_rms = [ 40, 115, 165, 215, 290, 350]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(168,180)
    file_in='L1544/L1544_NH3_11.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

    print("Now NH3(2,2)")
    a_rms = [  0,  20]
    b_rms = [ 10,  25]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(10,18)
    file_in='L1544/L1544_NH3_22.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

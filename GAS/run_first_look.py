from . import first_look
import numpy as np


def FirstLook_OrionA():
    print("Now NH3(1,1)")
    a_rms = [  0, 158, 315, 428, 530, 693, 751]
    b_rms = [ 60, 230, 327, 438, 604, 735, 760]
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
    a_rms = [  0, 115, 280, 385, 490, 655]
    b_rms = [ 80, 230, 345, 455, 625, 760]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(352,381)
    file_in='B18/B18_NH3_11.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now CCS")
    a_rms = [   0, 245]
    b_rms = [ 210, 490]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(225,243)
    file_in='B18/B18_CCS.fits'
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)
    #
    print("Now HC5N")
    a_rms = [  10, 245]
    b_rms = [ 210, 540]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(225,243)
    file_in='B18/B18_HC5N.fits'
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
    #
    print("Now NH3(2,2)")
    a_rms = [  0, 230, 460, 665]
    b_rms = [150, 380, 610, 820]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(380,460)
    file_in='Serpens_Aquila_test/Serpens_Aquila_NH3_22.fits'
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
    file_in='Serpens_Aquila_test/Serpens_Aquila_NH3_33.fits'
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
    file_in='Serpens_Aquila_test/Serpens_Aquila_C2S.fits'
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
    file_in='Serpens_Aquila_test/Serpens_Aquila_HC5N.fits'
    # 1st order polynomial
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
    a_rms = [  0, 260, 520, 730]
    b_rms = [150, 380, 610, 850]
    index_rms=first_look.create_index( a_rms, b_rms)
    index_peak=np.arange(380,520)
    file_in='NGC1333/NGC1333_NH3_22.fits'
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
    file_in='NGC1333/NGC1333_NH3_33.fits'
    # 1st order polynomial
    file_out=file_in.replace('.fits','_base1.fits')
    file_new=first_look.baseline( file_in, file_out, index_clean=index_rms, polyorder=1)
    first_look.peak_rms( file_new, index_rms=index_rms, index_peak=index_peak)

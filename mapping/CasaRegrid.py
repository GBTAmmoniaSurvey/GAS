filename='NGC1333_ABCDEFGH_NH3_11_all_blsub.fits'
importfits(fitsimage=filename,
           imagename='NGC1333_n11.im',
           overwrite=True)
importfits(fitsimage='NGC1333_ABCDEFGH_NH3_22_all_blsub.fits',
           imagename='NGC1333_n22.im',
           overwrite=True)
imregrid(imagename='NGC1333_n11.im',
         template='NGC1333_n22.im',
         output='NGC1333_n11.rg.im',
         axes=[0,1])
exportfits(imagename='NGC1333_n11.rg.im',
           fitsimage=filename.replace('blsub','blsub.rg'))




##########
First Look
##########

We have also created a first look pipeline that will allow us to baseline all the cubes and to create integrated intensity maps, peak temperature and first moment maps. 

This is carried out with wrapper functions for each observed region, which can be found here `First Look wrappers <../GAS.run_first_look>`, where the underlying first look functionality is defined in the `First Look subpackage <../GAS.first_look>`.

Go to the directory with the survey’s images and open an ipython session:
.. code-block:: python
    cd /lustre/pipeline/scratch/GAS/images/
    import GAS.run_first_look
    GAS.run_first_look.FirstLook_SerAqu(file_extension='_DR1')

Here the `file_extension` keyword is used to work with different versions of the data, e.g., DR1 and DR2. 

To know which first look functions are available just type:
.. code-block:: python
    import GAS.run_first_look
    GAS.run_first_look.<Press TAB>

*Important: The output from these functions are FITS files which are good enough to assess the quality of the data. A higher quality product is (will be) created using the results from the line fitting as input to determine an optimal mask.*

************************
What does it do for you?
************************

The First Look pipeline will take a FITS cube and remove the a polynomial baseline, for this you need to give the range of line free channels to be used and the polynomial order to be fit. Notice that we use a single channel range and polynomial order for the entire cube. Also, we use this same selection of channels to estimate he rms in each pixel. 

Also, the pipeline will determine what is the line peak brightness at each pixel, for this you need to provide a single range of channels where the emission is constrained (this is done to avoid noise spikes). We also calculate the *first look* integrated intensity map, using the same channel range as for the peak brightness. *Caution* this integrated intensity map only takes into account the channel range provided, which is OK for typical lines however it means that the satelite components will be missed from the NH\ :sub:`3` maps. The solution for this is provided by `update_NH3_moment0`.

All products (integrated intensity, rms, and peak brightness maps as well as the baselined cube) are exported as FITS files.

******************************
How to use it on non-GAS data?
******************************

If you are using this package for your own non-GAS data, then you need to do the following.

TODO

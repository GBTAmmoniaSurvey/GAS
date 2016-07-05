<<<<<<< HEAD
Line Fitting
============

The line fitting is done using ``pyspeckit``
=======
############
Line Fitting
############

The line fitting is done using `pyspeckit <http://pyspeckit.bitbucket.org>`_
We use the NH3 model included in ``pyspeckit`` to simultaneously fit NH\ :sub:`3` (1,1) and (2,2) to derive the parameters of NH3.
This includes the centroi velocity, velocity dispersion, excitation temperature, 
kinetic temperature, and column density.
Although we also observed NH\ :sub:`3` (3,3) we do not attempt to fit for the ortho-para ratio of NH\ :sub:`3`, by default we set the ortho-para ratio to 0, therefore, we only report the para-NH3 column density.


****************
Fitting GAS data
****************

We have setup a couple of convenience functions to do the line fitting. For this, it assumes that the files are stored in a directory named 'region_name'. 
For example, if we want to fit the data for OrionA, then we do the following:

  .. code-block:: python
  
    import GAS.PropertyMaps
    GAS.PropertyMaps.cubefit(region='OrionA', vmin=5.6, vmax=13.7, do_plot=False, snr_min=3.0, 
            multicore=40, file_extension='base_DR1', mask_function = None)

If you want to make some quick plots showing the fit results, then do the following

  .. code-block:: python
  
    import astropy.units as u
    import GAS.PropertyMaps
    GAS.PropertyMaps.plot_cubefit(region='OrionA', distance=450*u.pc, dvmin=0.05, dvmax=0.7, 
                 vcmin=5.7, vcmax=12.7, file_extension='base_DR1')

The result from this line fitting will most likely include pixels that don't have a good fit (because of the default value for the minimum snr), however, we only need to do this fit once and now we clean-up the parameter fit to keep only those that are reliable.
>>>>>>> GBTAmmoniaSurvey/master

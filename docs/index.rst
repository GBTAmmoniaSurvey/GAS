GAS documentation
=================

The ``GAS`` package provides a full data reduction pipeline of the data 
obtained by the GouldBelt Ammonia Survey (GAS) using the Green Bank Telescope. 

It provides a uniform framework to reduce each region observed in a consistent 
fashion. At the moment, the GAS pipeline does the following:

- It runs the official calibration pipeline from GBT 
  (https://github.com/nrao/gbt-pipeline/) on the raw data. This is done through 
  an wrapper.
- Each spectral window and region (as defined in the observing log) are stored 
  in a different directory.
- The imaging is done with our gridder, which generates a data cube using all 
  the observations for a given window.
- A first look functionality. This allows for a quick look of the data, useful 
  for quality assesment.
- Line fitting for all data is done using ``pyspeckit``. In the case of NH3 
  the proper hyperfine model is used, while for the other lines a single 
  Gaussian model is used.


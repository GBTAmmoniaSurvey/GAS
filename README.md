About
-----
-----

This repository includes the scripts used to prepare the GBT Ammonia Survey (GAS, PIs: Jaime E Pineda and Rachel Friesen). 

It also facilitates the data reduction and analysis of all the data produced by the survey. 
We have developed a pipeline that automatically reduced the survey data, with calibration done 
using the gbt-pipeline, imaging carried out with a custom made gridding package, and line-fitting 
carried out using pyspeckit.

More information is available in the documentation, avaliable [online at readthedocs.io](http://gas.readthedocs.io).

Credits
-------
-------

This is developed by:
* Jaime E Pineda ([@jpinedaf](http://github.com/jpinedaf))
* Erik Rosolowsky ([@low-sky](http://github.com/low-sky))
* Rachel Friesen ([@rfriesen](http://github.com/rfriesen))
* Adam Ginsburg ([@keflavich](http://github.com/keflavich))

Dependencies
------------
------------

* astropy (>=1.1)
* pyspeckit (>0.1.18.1)
* matplotlib (>=1.5.1)
* aplpy (>=1.0)
* spectral_cube (>=0.3.0)
* radio_beam
* scikit-image (>=0.12.3)

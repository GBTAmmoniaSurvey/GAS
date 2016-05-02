.. _section_linefit:

Line Fitting
============

The line fitting is done using `pyspeckit <htpp://pyspeckit.bitbucket.org>`_
We use the NH3 model included in ``pyspeckit`` to simultaneously fit NH3 (1,1) and (2,2) to derive the parameters of NH3.
This includes the centroi velocity, velocity dispersion, excitation temperature, 
kinetic temperature, and column density.
Although we also observed NH3(3,3) we do not attempt to fit for the ortho-para ratio of NH3, by default we set the ortho-para ratio to 0, therefore, we only report the para-NH3 column density.



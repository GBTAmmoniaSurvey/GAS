# Licensed under a MIT style license - see LICENSE

"""
This is a package to process the GBT Ammonia Survey (GAS) data. It includes a calibration funtion 
that (wrapper of gbtpipeline), an imaging function which creates the reduced data cubes, and the 
first look analysis tools. It is setup to run at Green Bank. This package relies in Astropy.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:
    import gasPipeline
    import gridregion
    import run_grid_regions
    import first_look
    import run_first_look
    import PropertyMaps
    import gasBinning
    import voronoi_2d_binning
    pass

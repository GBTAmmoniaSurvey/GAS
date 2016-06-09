###########
Calibration
###########

The calibration is done in a region by region basis, as listed in the observations summary. Go to the directory of the region and open an ipython session,

  .. code-block:: python
  
    cd /lustre/pipeline/scratch/GAS/region_name
    import GAS
    GAS.gasPipeline.wrapper( region=region_name, startdate=’2015-11-10’, enddate=’2015-12-10’)  

if redoing this step and need to download the latest observation summary file, remove the .csv file in the directory you are running it first).

The first command will load the calibration pipeline functions, while the second downloads an up-to-date version of the observations summary and reduces all spectral windows and beams observed for the requested region. 
Also, there is a function that runs the data reduction for all regions, just run the following commands,

  .. code-block:: python
  
    cd /lustre/pipeline/scratch/GAS/
    import GAS
    GAS.gasPipeline.reduceAll()

The pipeline puts the calibrated data in named directories for each region (i.e., for Orion A, the calibrated data are in /lustre/pipeline/scratch/GAS/OrionA). Within the named directories, there will be folders for each molecular line in the observing setup. 

*Important: Check against the observing logs to see if all the sessions/scans are being reduced for all lines. If not, take note and then we can check what to do, or if we need to reduce them manually.*
#######
Imaging
#######

Go to the directory with the survey’s images for the region to generate::

    >>> cd /lustre/pipeline/scratch/GAS/images/region_name

and run the function to image the region want. This is an example for Serpens_Aquila::

  .. code-block:: python
  
    cd /lustre/pipeline/scratch/GAS/images/Serpens_Aquila
    import GAS
    GAS.run_grid_regions.grid_SerAqu()

The pipeline puts the cubes in an ‘images’ directory, with subdirectories for each region (i.e., for Orion A, the calibrated data are in /lustre/pipeline/scratch/GAS/images/OrionA). 

Please do check what is the proper channel range needed to include all the data, while avoiding too many empty channels. Once you know the prefered range then modify the grid function in your local distribution of the pipeline and commit. Don’t modify the files at Green Bank computers.

To know which imaging functions are available just type

  .. code-block:: python
  
    import GAS
    GAS.run_grid_regions.<Press TAB>
    GAS.run_grid_regions.fits          GAS.run_grid_regions.grid_L1688    GAS.run_grid_regions.grid_OrionB
    GAS.run_grid_regions.grid_B18      GAS.run_grid_regions.grid_NGC1333  GAS.run_grid_regions.grid_SerAqu
    GAS.run_grid_regions.grid_L1455    GAS.run_grid_regions.grid_OrionA   GAS.run_grid_regions.gridregion


First Look
==========

We have also created a first look pipeline that will allow us to baseline all the cubes and to create integrated intensity maps, peak temperature and first moment maps. 

This is carried out with wrapper functions for each observed region, which can be found here `First Look wrappers <../GAS.run_first_look>`, where the underlying first look functionality is defined in the `First Look subpackage <../GAS.first_look>`.

Go to the directory with the survey’s images and open an ipython session:
    >>> cd /lustre/pipeline/scratch/GAS/images/
    >>> ipython
    >>> import GAS.run_first_look
    >>> GAS.run_first_look.FirstLook_SerAqu(file_extension='_DR1')

Here the `file_extension` keyword is used to work with different versions of the data, e.g., DR1 and DR2. 

To know which first look functions are available just type:
    >>> ipython
    >>> import GAS.run_first_look
    >>> GAS.run_first_look.<Press TAB>

*Important: The output from these functions are FITS files which are good enough to assess the quality of the data. A higher quality product is (will be) created using the results from the line fitting as input to determine an optimal mask.*


How to use it on non-GAS data?
------------------------------

If you are using this package for your own non-GAS data, then you need to do the following.

TODO
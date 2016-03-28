First Look
==========

We have also created a first look pipeline that will allow us to baseline all the cubes and to create integrated intensity maps, peak temperature and first moment maps. 

This is carried out with wrapper functions for each observed region, which can be found here :doc:`First Look wrappers <GAS.run_first_look>`, where the underlying first look functionality is defined in the :doc:`First Look subpackage <GAS.first_look>`.

Go to the directory with the survey’s images and open an ipython session
code-block::
    >>> cd /lustre/pipeline/scratch/GAS/images/
    >>> ipython
    >>> import GAS
    >>> GAS.run_first_look.FirstLook_SerAqu()

To know which first look functions are available just type: ::
code-block::
    >>> ipython
    >>> import GAS
    >>> GAS.run_first_look.<Press TAB>

*Important: The output from these functions are FITS files which are good enough to assess the quality of the data. A higher quality product is (will be) created using the results from the line fitting as input to determine an optimal mask.*

import gridregion
import os

cwd = os.getcwd()
test_dir='grid_data_test'
test_out='test_grid_output'
os.chdir( test_out)
gridregion.griddata( rootdir = cwd, region = '', dirname = test_dir)

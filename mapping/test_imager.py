import gridregion
import os

cwd = os.getcwd()
test_dir='grid_data_test'

os.chdir( test_dir)
gridregion.griddata( rootdir = cwd, region = '', dirname = test_dir)

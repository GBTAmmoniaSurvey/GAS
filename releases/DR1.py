dr_id = 'DR1'

def install_and_import(package, url):
    # http://stackoverflow.com/questions/12332975/installing-python-module-within-code
    import importlib
    try:
        importlib.import_module(package)
    except ImportError:
        import pip
        pip.main(['install', package, '-t', 'GAS_{dr}'.format(dr=dr_id)])
    finally:
        globals()[package] = importlib.import_module(package)


install_and_import('GAS', 'git+https://github.com/GBTAmmoniaSurvey/GAS.git@c9b5ab260b1df81f823edbe64de4230a0f67dd09#egg=GAS')
#install_and_import('GAS', 'git+https://github.com/GBTAmmoniaSurvey/GAS.git@GAS_{dr}#egg=GAS'.format(dr=dr_id))

# Check that the appropriate version has been installed
assert GAS.__version__ == dr_id, "Wrong GAS version installed: {0} instead of {1}".format(GAS.__version__, dr_id)

release_version = 'DR1'
file_extension = '_{0}'.format(release_version)
#
if GAS.run_grid_region.grid_B1(release=release_version):
    GAS.run_first_look.FirstLook_B1(file_extension=file_extension)
#
if GAS.run_grid_region.grid_B1E(release=release_version):
    GAS.run_first_look.FirstLook_B1E(file_extension=file_extension)
#
if GAS.run_grid_region.grid_B18(release=release_version):
    GAS.run_first_look.FirstLook_B18(file_extension=file_extension)
#
if GAS.run_grid_region.grid_B59(release=release_version):
    GAS.run_first_look.FirstLook_B59(file_extension=file_extension)
#
if GAS.run_grid_region.grid_Cepheus(release=release_version):
    GAS.run_first_look.FirstLook_Cepheus(file_extension=file_extension)
#
if GAS.run_grid_region.grid_IC348(release=release_version):
    GAS.run_first_look.FirstLook_IC348(file_extension=file_extension)
#
if GAS.run_grid_region.grid_IC5146(release=release_version):
    GAS.run_first_look.FirstLook_IC5146(file_extension=file_extension)
#
if GAS.run_grid_region.grid_L1451(release=release_version):
    GAS.run_first_look.FirstLook_L1451(file_extension=file_extension)
#
if GAS.run_grid_region.grid_L1455(release=release_version):
    GAS.run_first_look.FirstLook_L1455(file_extension=file_extension)
#
if GAS.run_grid_region.grid_L1688(release=release_version):
    GAS.run_first_look.FirstLook_L1688(file_extension=file_extension)
#
if GAS.run_grid_region.grid_L1689(release=release_version):
    GAS.run_first_look.FirstLook_L1689(file_extension=file_extension)
#
if GAS.run_grid_region.grid_NGC1333(release=release_version):
    GAS.run_first_look.FirstLook_NGC1333(file_extension=file_extension)
#
if GAS.run_grid_region.grid_OrionA(release=release_version):
    GAS.run_first_look.FirstLook_OrionA(file_extension=file_extension)
#
if GAS.run_grid_region.grid_OrionB_NGC2023_2024(release=release_version):
    GAS.run_first_look.FirstLook_OrionB_NGC2023_2024(file_extension=file_extension)
#
if GAS.run_grid_region.grid_OrionB_NGC2068_2071(release=release_version):
    GAS.run_first_look.FirstLook_OrionB_NGC2068_2071(file_extension=file_extension)
#
if GAS.run_grid_region.grid_SerAqu(release=release_version):
    GAS.run_first_look.FirstLook_SerAqu(file_extension=file_extension)
#
if GAS.run_grid_region.grid_TMC1(release=release_version):
    GAS.run_first_look.FirstLook_TMC1(file_extension=file_extension)
#
if GAS.run_grid_region.grid_HC2(release=release_version):
    GAS.run_first_look.FirstLook_HC2(file_extension=file_extension)

dr_id = 'DR1'

def install_and_import(package, url):
    # http://stackoverflow.com/questions/12332975/installing-python-module-within-code
    import importlib
    try:
        importlib.import_module(package)
    except ImportError:
        import pip
        pip.main(['install', '-e', package, '-t', 'GAS_{dr}'.format(dr=dr_id)])
    finally:
        globals()[package] = importlib.import_module(package)


install_and_import('GAS', 'git+https://github.com/GBTAmmoniaSurvey/GAS.git@c9b5ab260b1df81f823edbe64de4230a0f67dd09#egg=GAS')
#install_and_import('GAS', 'git+https://github.com/GBTAmmoniaSurvey/GAS.git@GAS_{dr}#egg=GAS'.format(dr=dr_id))

GAS.run_grid_region.grid_B1()
GAS.run_grid_region.grid_B18()
GAS.run_grid_region.grid_B59()
GAS.run_grid_region.grid_Cepheus()
GAS.run_grid_region.grid_IC348()
GAS.run_grid_region.grid_IC5146()
GAS.run_grid_region.grid_L1451()
GAS.run_grid_region.grid_L1455()
GAS.run_grid_region.grid_L1688()
GAS.run_grid_region.grid_NGC1333()
GAS.run_grid_region.grid_OrionA()
GAS.run_grid_region.grid_OrionB()
GAS.run_grid_region.grid_SerAqu()

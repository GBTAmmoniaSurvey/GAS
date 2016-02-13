from __future__ import absolute_import

import os

def get_package_data():
    paths = [os.path.join('data','ObservationLog_DR1.csv'),
             os.path.join('data','RegionCatalog_DR1.csv')]
    return {'GAS': paths}

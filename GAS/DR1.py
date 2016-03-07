from . import gasPipeline
from . import run_first_look
from . import PropertyMaps

def DR1():
    gasPipeline.reduceAll(release='DR1')
    run_first_look.FirstLook(release='DR1')


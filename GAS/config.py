import astropy.units as u
import numpy as np

L1688 = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
         'distance' : 137.3*u.pc,
         'sigma_mult' : 1.
         }

B18 = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
       'distance' : 135*u.pc,
       'sigma_mult' : 1.
       }

NGC1333 = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
           'distance' : 250.*u.pc,
           'sigma_mult' : 1.
         }

OrionA = {'scalebar_size': 0.5*u.pc, 'scalebar_pos': 'bottom right',
          'distance' : 414.*u.pc,
          'sigma_mult' : 1.5
         }
#Gobelins distance
OrionA_S = {'scalebar_size': 0.5*u.pc, 'scalebar_pos': 'bottom right',
            'distance' : 428.*u.pc, 
            'sigma_mult' : 1.
         }

B1 = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
      'distance' : 250.*u.pc,
      'sigma_mult' : 1.
         }

CrAeast = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
           'distance' : 130.*u.pc,
           'sigma_mult' : 1.
         }

HC2 = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
       'distance' : 135.*u.pc,
       'sigma_mult' : 1.
         }

B1E = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
       'distance' : 250.*u.pc,
       'sigma_mult' : 1.
       }

CrAwest = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
           'distance' : 130.*u.pc,
           'sigma_mult' : 1.}

B59 = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
       'distance' : 140.*u.pc,
       'sigma_mult' : 1.}

Cepheus_L1228 = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
                 'distance' : 300.*u.pc,
                 'sigma_mult' : 1.}

Cepheus_L1251 = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
                 'distance' : 300.*u.pc,
                 'sigma_mult' : 1.}

IC348 = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
         'distance' : 316.*u.pc,
         'sigma_mult' : 1.}

IC5146 = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
          'distance' : 460.*u.pc,
          'sigma_mult' : 1.}

OrionB_NGC2023_2024 = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
                       'distance' : 420.*u.pc,
                       'sigma_mult' : 1.}

OrionB_NGC2068_2071 = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
                       'distance' : 388.*u.pc,
                       'sigma_mult' : 1.}

Perseus = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
           'distance' : 250.*u.pc,
           'sigma_mult' : 1.}

Pipe_Core40 = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
               'distance' : 450.*u.pc,
               'sigma_mult' : 1.}

Serpens_MWC297 = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
                  'distance': 437.*u.pc,
                  'sigma_mult' : 1.}

Serpens_Aquila = {'scalebar_size': 0.1*u.pc, 'scalebar_pos': 'bottom right',
                  'distance': 437.*u.pc,
                  'sigma_mult' : 1.}

plottingDictionary = {"L1688" : L1688, "B18" : B18, "NGC1333" : NGC1333, 
                      "OrionA" : OrionA, "OrionA_S" : OrionA_S, "B1" : B1, 
                      "HC2" : HC2, "Pipe_Core40" : Pipe_Core40, "Perseus" : Perseus,
                      "OrionB_NGC2068-2071" : OrionB_NGC2068_2071, 
                      "OrionB_NGC2023-2024" : OrionB_NGC2023_2024,
                      "IC5146" : IC5146, "IC348" : IC348, "Cepheus_L1228" : Cepheus_L1228, 
                      "Cepheus_L1251" : Cepheus_L1251, "B59" : B59, "CrAeast" : CrAeast,
                      "CrAwest" : CrAwest, "B1E" : B1E, "Serpens_Aquila" : Serpens_Aquila,
                      "Serpens_MWC297" : Serpens_MWC297}

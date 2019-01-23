# -*- coding: utf-8 -*-

#child class of SpecData. Initializes SpecData to create the amplitudes and
# axes.

# this is where we will add the Mims ENDOR methods like ln, baseline, exp
# deglitch, and flip to the initialization.

# we can add a method for the calculation of ENDOR freq x.
from Spec_Data import Spec_Data

##############################################################################
##############################################################################

class MimsENDOR(Spec_Data):


    def __init__(self, metadata, full_file_name):
        super().__init__( metadata, full_file_name)

        #this knows self.phased from SpecData. It is the phased data that is ready
        # to be processed further.
        print(len(self.phased))
        # this is just to check that it knows what self.amp is.


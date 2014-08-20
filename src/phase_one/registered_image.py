"""Class that registers two images and stores info about the registration.
"""

import os
import subprocess

class RegisteredImage:

    def __init__(self, reference, movable):

        # TODO validation
        self._reference = reference
        self._movable = movable

        self._robust_reg_options = ['--satit',
                                    '--iscale']

    def register(self):

        cmd = ['mri_robust_register',
               '--mov', self._movable,
               '--dst', self._reference,
               '--lta', os.path.splitext(self._movable)[0] + '.lta',
               '--mapmov', os.path.splitext(self._movable)[0] + '_map.nii',
               '--weights', os.path.splitext(self._movable)[0] + '_weights.nii',
              ] + self._robust_reg_options

        try:
            out = subprocess.check_output(cmd)
        except:
            print "Failure executing command:"
            print cmd
            return False

        return True

"""Class that registers two images and stores info about the registration.
"""

import os
import sys
import subprocess

class RegisteredImage:

    _reg_prog = 'mri_robust_register'

    @classmethod
    def check_environment(cls):
        cmd = [cls._reg_prog]

        try:
            check_proc = subprocess.Popen(cmd,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE)
            out, err = check_proc.communicate()

            if (check_proc.returncode == 1 and
                out.find(cls._reg_prog) != -1):
                return True
            else:
                print "Unexpected output while executing %s" % cls._reg_prog
                return False

        except(OSError):
            print "Can't find %s, make sure it's in the $PATH?" % cls._reg_prog
            return False


    def __init__(self, reference, movable):
        self._reference = reference
        self._movable = movable

        self._robust_reg_options = ['--satit',
                                    '--iscale']

    def register(self):

        out_stem = self._movable.split('.')[0]
        cmd = [self._reg_prog,
               '--mov', self._movable,
               '--dst', self._reference,
               '--lta', out_stem + '.lta',
               '--mapmov', out_stem + '_map.nii.gz',
               '--weights', out_stem + '_weights.nii.gz',
              ] + self._robust_reg_options

        try:
            reg_proc = subprocess.Popen(cmd,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
            out, err = reg_proc.communicate()

            if reg_proc.returncode == 0:
                with open(out_stem + '.log', 'w') as log:
                    log.write(err)
            else:
                print "Failure executing command:"
                print cmd
                print err
                return False

        except(OSError):
            print ("Error executing subprocess for registration, perhaps %s is "
                   "not in the $PATH?" % self._reg_prog)
            return False

        return True

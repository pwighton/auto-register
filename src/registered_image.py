#!/usr/bin/env python

"""Class that registers two images and stores info about the registration.
"""

import argparse
import nibabel as nb
import numpy as np
import os
import sys
import subprocess
import re

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

    @classmethod
    def read_transform_file(cls, filename):
        with open(filename) as f:
            transform_line_re = re.compile('([0-9.e\-\+]+\s){4}')
            transform = ''
            num_consecutive_transform_lines = 0
            for line in f:
                if transform_line_re.match(line):
                    transform += line.strip() + ' '
                    num_consecutive_transform_lines += 1
                else:
                    if num_consecutive_transform_lines == 4:
                        break
                    else:
                        num_consecutive_transform_lines = 0
                        transform = ''

            # validate
            if num_consecutive_transform_lines != 4:
                print "Transform file ended unexpectedly"
                return ''

            if ([float(s) for s in transform.split()[-4:]] !=
                [0.0, 0.0, 0.0, 1.0]):
                print "Transform file parse error"
                print [float(s) for s in transform.split()[-4:]]

                return ''

        return ' '.join(transform.split()[:-4])

    def __init__(self, reference, movable):
        self._reference = reference
        self._movable = movable

        self._robust_reg_options = ['--satit',
                                    '--iscale']

        self._transform_file = None

    def get_transform_filename(self):
        return self._transform_file

    def get_transform(self):
        return RegisteredImage.read_transform_file(self._transform_file)

    def register(self):

        out_stem = '.'.join(self._movable.split('.')[:-1])

        # also handle .nii.gz files
        if out_stem.endswith('.nii'):
            out_stem = out_stem[:-4]

        self._transform_file = out_stem + '.lta'

        cmd = [self._reg_prog,
               '--mov', self._movable,
               '--dst', self._reference,
               '--lta', self._transform_file,
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


def main(argv):
    def verifyPathExists(path):
        if not os.path.exists(path):
            raise ValueError("%s does not exist" % path)
        return path

    def verifyDirExists(path):
        if not os.path.isdir(path):
            raise ValueError("%s does not exist" % path)
        return path

    # parse args
    parser = argparse.ArgumentParser()
    parser.add_argument('reference', type=verifyPathExists,
                        help='Path to nifti image to use as the reference '
                        'image for registration')
    parser.add_argument('movable', type=verifyPathExists,
                        help='Path to nifti image to use as the movable '
                        'image for registration')

    args = parser.parse_args()

    ri = RegisteredImage(args.reference, args.movable)

    print "Registering %s to %s" % (args.movable, args.reference)

    if ri.register():
        print "Registration completed"
    else:
        print "Registration failed"


if __name__ == "__main__":
    sys.exit(main(sys.argv))

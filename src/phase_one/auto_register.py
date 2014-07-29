#!/usr/bin/env python

import argparse
import os
import sys

import ext.receive_nii

class AutoRegister:

    def __init__(self):
        self._should_shutdown = False
        self._receiver = receive_nii.ImageReceiver(args)
        self._receiver.start()

    def run(self):
        while self._receiver.is_running() and not self._should_shutdown:
            filename = self._receiver.get_next_filename()
            if filename is not None:
                print "saw file %s" % filename

                # TODO handle the file


def main(args):
    AutoRegister ar
    ar.run()
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))

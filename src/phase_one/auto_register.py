#!/usr/bin/env python

import argparse
import os
import sys
import thread
import traceback

from ext import receive_nii
from util import TerminalInput

class AutoRegister:

    def __init__(self, args):
        self._should_shutdown = False

        # unused but required reciever args
        args.four_dimensional = False
        args.single_series = False

        self._receiver = receive_nii.ImageReceiver(args)
        self._term_input = TerminalInput()

        self._receiver.start()
        self._term_input.start()

    def check_for_input(self):
        c = self._term_input.get_char()
        if c == None:
            return

        if c == 'q':
            self.shutdown()


    def run(self):
        print "Running AutoRegister, 'q' to quit"

        while self._receiver.is_running() and not self._should_shutdown:
            self.check_for_input()

            filename = self._receiver.get_next_filename()
            if filename is not None:
                print "saw file %s" % filename

                # TODO handle the file

    def shutdown(self):
        self._should_shutdown = True
        self._term_input.stop()
        self._receiver.stop()
        print "Shuting down"

def main(args):
    def verifyPathExists(path):
        if not os.path.exists(path):
            raise ArgumentError()
        return path

    # parse args
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--save_directory', type=verifyPathExists,
                        help='Directory in which to save images and data')
    parser.add_argument('-H', '--host', default='localhost',
                        help='Address of the scanner from which to listen '
                        'for images [localhost]')
    parser.add_argument('-p', '--port', default=15001, type=int,
                        help='On the scanner address from which to listen '
                        'for images [15001]')

    ar = AutoRegister(parser.parse_args())
    ar.run()
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))

#!/usr/bin/env python

"""Main file and class for the autoregister application.
"""

import argparse
import os
import sys
import traceback

from image_receiver import ImageReceiver
from registered_image import RegisteredImage
from transform_sender import TransformSender
from util import TerminalInput

class AutoRegister(object):

    def __init__(self, args):
        """Initialize the autoregister application and helper modules.
        """
        # validate environment
        if not RegisteredImage.check_environment():
            raise ValueError("Environment check failed")

        # validate args
        if args.reference is None and not args.first:
            raise ValueError("One of --reference or --first must be set")
        elif args.reference is not None and args.first:
            raise ValueError("Both --reference and --first cannot be set")

        self._reference = args.reference
        self._should_shutdown = False

        # unused but required reciever args
        args.four_dimensional = False
        args.single_series = False
        self._image_receiver = ImageReceiver(args)

        self._transform_sender = TransformSender(15001)

        self._term_input = TerminalInput(disabled=args.no_terminal)

    def check_for_input(self):
        """Return the last character input, or None. If 'q' is seen, the
        autoregister application shuts down.
        """
        c = self._term_input.get_char()
        if c == None:
            return None

        if c == 'q':
            self.shutdown()
            return None

    def run(self):
        """Main loop of the autoregister application.
        """
        self._image_receiver.start()
        self._transform_sender.start()
        self._term_input.start()

        print "Running AutoRegister, 'q' to quit"

        while not self._should_shutdown:

            filename = self._image_receiver.get_next_filename()
            if filename is not None:
                if self._reference is None: # need a reference
                    self._reference = filename
                    print "Using reference: %s" % filename
                else: # register
                    reg_image = RegisteredImage(self._reference, filename)
                    print "Registering %s to the reference" % filename
                    if reg_image.register():
                        print "Registration complete"

                        if self._transform_sender.send(
                                reg_image.get_transform_filename()):
                            print "Transform ready to send"
                        else:
                            print "Failed to prepare transform for sending"

            # must be the last task in the mainloop to handle shutdown
            # properly
            self.check_for_input()


    def shutdown(self):
        """Shutdown the autoregister application. Stops the mainloop and
        tears down helper modules.
        """
        print "Shuting down"
        self._should_shutdown = True
        self._term_input.stop()
        self._image_receiver.stop()
        self._transform_sender.stop()

def main(args):
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
    parser.add_argument('-s', '--save_directory', type=verifyDirExists,
                        required=True,
                        help='Directory in which to save images and data')
    parser.add_argument('-r', '--reference', type=verifyPathExists,
                        help='Path to nifti image to use as the reference '
                        'image for registration')
    parser.add_argument('-f', '--first', action="store_true",
                        help='Take first image received as the reference '
                        'image for registration')
    parser.add_argument('-H', '--host', default='localhost',
                        help='Address of the scanner from which to listen '
                        'for images [localhost]')
    parser.add_argument('-p', '--port', default=15000, type=int,
                        help='On the scanner address from which to listen '
                        'for images [15000]')
    parser.add_argument('-T', '--no-terminal', action='store_true',
                        help='Do not listen for terminal input (helpful '
                        'for debugging)')

    ar = AutoRegister(parser.parse_args())
    ar.run()
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))

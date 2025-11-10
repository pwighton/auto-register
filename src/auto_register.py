#!/usr/bin/env python2

"""Main file and class for the autoregister application.
"""

import argparse
import os
import sys
import traceback
import numpy as np
import signal
import shutil

from image_manager import ImageManager
from registered_image import RegisteredImage
from transform_sender import TransformSender
from terminal_input import TerminalInput
from np4x4_utils import *

class AutoRegister(object):

    def __init__(self, args):
        """Initialize the autoregister application and helper modules.
        """
        # Set up signal handler for graceful shutdown
        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)
        
        # validate environment
        if not RegisteredImage.check_environment():
            raise ValueError("RegisteredImage Environment check failed")

        # validate args
        if args.reference is None and not args.first and args.transform is None:
            raise ValueError("One of --reference, --first or --transform must be set")
        elif args.reference is not None and args.first:
            raise ValueError("Both --reference and --first cannot be set")

        self._reference = args.reference
        if self._reference is not None:
            print "Using reference: %s" % self._reference

        self._should_shutdown = False
        self._disable_default_areg = getattr(args, 'disable_default_areg', False)

        self._image_manager = ImageManager(args)

        # If the image manager is looking for dicoms in a dir, dcm2niix needs to be
        # in PATH
        if not self._image_manager.check_environment():
          raise ValueError("ImageManager Environment check failed")

        self._transform_sender = TransformSender(args.host, 15001, args.transform)
        self._term_input = TerminalInput(disabled=args.no_terminal)

        self._prescription_transform = string_to_np4x4(args.prescription)
        if self._prescription_transform is not None:
            print "Prescription matrix has been defined.  All registrations will be multiplied by"
            print self._prescription_transform
        self._last_transform = None
        self._resend_mode = args.resend
        if args.transform is not None:
            self._resend_mode = True

    def _signal_handler(self, signum, frame):
        """Handle interrupt signals (Ctrl+C) gracefully."""
        print "\nReceived interrupt signal, shutting down gracefully..."
        self.shutdown()
        
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
        self._image_manager.start()
        self._transform_sender.start()
        self._term_input.start()

        print "Running AutoRegister, 'q' to quit"
        if len(self._transform_sender._transforms_to_send) > 0:
          print "Will manually send the following transform: "
          print "  ", self._transform_sender._transforms_to_send[0]

        while not self._should_shutdown:

            self._transform_sender.set_state("checking")
            filename = self._image_manager.get_next_filename()
            if filename is not None:
                if self._disable_default_areg:
                    # Default autoregistration is disabled, just report the file
                    print "Received file: %s (default autoregistration disabled)" % filename
                elif self._reference is None: # need a reference
                    self._reference = filename
                    print "Using reference: %s" % filename
                else: # register
                    self._transform_sender.set_state("registering")

                    reg_image = RegisteredImage(self._reference, filename,
                                                verbose=True)

                    def send_last_transform():
                        if self._transform_sender.send(self._last_transform):
                            print "Transform ready to send"
                        else:
                            print "Failed to prepare transform for sending"

                    if self._resend_mode and self._last_transform is not None:
                        print ("WARNING: in resend mode, just resending the "
                               "transform we already have")
                        send_last_transform()
                    else:
                        print "Registering %s to the reference" % filename
                        if reg_image.register():
                            print "Registration complete"
                            print "Registration Transform:"
                            print reg_image.get_transform()

                            # If --prescription was defined, multiply the registration with the prescription
                            if self._prescription_transform is None:
                                self._last_transform = reg_image.get_transform()
                            else:
                                print "Prescription Transform:"
                                print self._prescription_transform
                                reg_img_transform = string_to_np4x4(reg_image.get_transform())
                                last_transform = np.matmul(reg_img_transform, self._prescription_transform)
                                self._last_transform = np4x4_to_string(last_transform)
                            print "Transform to send to scanner:"
                            print string_to_np4x4(self._last_transform)
                            send_last_transform()

            self._transform_sender.clear_state()

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
        self._image_manager.stop()
        self._transform_sender.stop()

def main(args):
    """Main entry point
    """

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
    parser.add_argument('-R', '--resend', action="store_true",
                        help='After a registration is computed, just resend '
                        'that registration each time a new image is received '
                        '(useful for minimizing the time impact of persistent '
                        'timeout-based image recon errors)')
    parser.add_argument('--input-mode', choices=['socket', 'directory'], 
                        default='socket',
                        help='Input mode: socket (TCP) or directory (monitor for DICOM files) [socket]')
    parser.add_argument('--watch-directory', type=verifyDirExists,
                        help='Directory to monitor for DICOM files (required when input-mode=directory)')
    parser.add_argument('--poll-interval', type=float, default=2.0,
                        help='Seconds between directory scans when in directory mode [2.0]')
    parser.add_argument('--stabilize-wait', type=float, default=2.0,
                        help='Seconds to wait for file size to stabilize before processing when in directory mode [2.0]')
    parser.add_argument('--dicom-filter', type=verifyPathExists,
                        help='JSON file containing DICOM tag/value pairs to filter incoming files (directory mode only)')
    parser.add_argument('--command', action='append', dest='commands',
                        help='Arbitrary command to execute when a new DICOM is detected (can be specified multiple times)')
    parser.add_argument('--command-filter', action='append', type=verifyPathExists, dest='command_filters',
                        help='JSON file containing DICOM tag/value pairs to filter for corresponding --command (can be specified multiple times)')
    parser.add_argument('--disable-default-areg', action='store_true',
                        help='Disable the default behavior of calling mri_robust_register')
    parser.add_argument('-H', '--host', default='0.0.0.0',
                        help='Address of the scanner from which to listen '
                        'for images [0.0.0.0]')
    parser.add_argument('-p', '--port', default=15000, type=int,
                        help='On the scanner address from which to listen '
                        'for images [15000]')
    parser.add_argument('-T', '--no-terminal', action='store_true',
                        help='Do not listen for terminal input (helpful '
                        'for debugging)')
    parser.add_argument('-syn', '--synthstrip', action='store_true',
                        default=False,
                        help='Run mri_synthstrip on incomming niftis (mri_synthstrip must be in $PATH)')
    parser.add_argument('-trans', '--transform', type=str,
                        help='Manually specify a transformation matrix to send to scanner (string with 16 floats)',
                        default=None)
    parser.add_argument('-prescrip', '--prescription', type=str,
                        help='Specify a prescription matrix (text file with 16 floats; LPS) Registration matrix will be multiplied by this matrix (M_regsiter * M_prescription) and the result will be sent to the scanner')
    parser.add_argument('--expect-preheader', action='store_true', default=True,
                        help='Expect vsend to send an 8-byte preheader (default: True)')
    parser.add_argument('--no-expect-preheader', action='store_false', dest='expect_preheader',
                        help='Do not expect vsend to send an 8-byte preheader when in socket mode')
    
    args = parser.parse_args()
    print "Command line args: ", args

    # Validate incompatible flag combinations
    if args.disable_default_areg and args.dicom_filter is not None:
        raise ValueError("--disable-default-areg and --dicom-filter cannot be used together")

    # Pair up commands with filters
    command_filter_pairs = []
    if args.commands is not None:
        num_commands = len(args.commands)
        num_filters = len(args.command_filters) if args.command_filters is not None else 0

        for i in range(num_commands):
            command = args.commands[i]
            filter_path = args.command_filters[i] if i < num_filters else None
            command_filter_pairs.append((command, filter_path))

        # Validate that all commands exist
        for command, _ in command_filter_pairs:
            # Extract the first token (the actual command/executable)
            cmd_parts = command.split()
            if len(cmd_parts) == 0:
                raise ValueError("Empty command specified")

            executable = cmd_parts[0]

            # Check if it's an absolute path
            if os.path.isabs(executable):
                if not os.path.exists(executable):
                    raise ValueError("Command not found: %s" % executable)
            else:
                # Check in PATH
                if shutil.which(executable) is None:
                    raise ValueError("Command not found in PATH: %s" % executable)

    # Store the paired commands in args for later use
    args.command_filter_pairs = command_filter_pairs

    try:
      ar = AutoRegister(args)
      ar.run()
      return 0
    except KeyboardInterrupt:
        print "\nInterrupted by user"
        return 1
    except Exception as e:
        print "Error: %s" % str(e)
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main(sys.argv))

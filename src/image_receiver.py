#!/usr/bin/env python

"""Receive images sent to a TCP server. Also save the images and make
the filename available through an external interface.

"""

import argparse
from collections import namedtuple
import glob
import os
import sys
from time import sleep
import threading

from tcpip_server import ThreadedTCPServer

import external_image
import nibabel as nb
import numpy as np


class ImageReceiver(object):
    """Run a TCP server, receive images, and save them.
    """

    def __init__(self, args):
        """Store arguments, determine first available image name so we don't
        overwrite existing images, and create a template image header for saving.
        """

        self.host = args.host
        self.port = args.port
        self.server = None

        self.save_location = args.save_directory
        self.save_volume_index = 0

        # find existing volumes so we don't overwrite them
        existing = glob.glob(os.path.join(
            self.save_location, "img-[0-9][0-9][0-9][0-9][0-9].nii.gz"))
        self.save_volume_index = len(existing)

        self.ei = external_image.ExternalImage("ExternalImageHeader")

        self.mutex = threading.Lock()
        self.filename_stack = []

    def stop(self):
        """Stop listening for incoming images.
        """
        if self.server is not None:
            self.server.shutdown()
            self.server = None

        print "Image receiver stopped"

    def start(self):
        """Start the server to listen for incoming images.
        """
        if self.server is not None and self.server.is_running():
            raise RuntimeError('Server already running')

        server = ThreadedTCPServer((self.host, self.port), self.process_data)
        ip, port = server.server_address
        print "Image receiver running at %s on port %d" % (ip, port)
        self.server = server

    def get_next_filename(self):
        """If there is a new image file available, return it and remove the
        filename from those available.
        """

        filename = None
        self.mutex.acquire()
        if len(self.filename_stack) > 0:
            filename = self.filename_stack.pop()
        self.mutex.release()

        return filename

    def is_running(self):
        """Get whether the server is running.
        """

        return self.server.is_running()

    def process_data(self, sock):
        """Callback to receive image data when it arrives.
        """

        in_bytes = sock.recv(self.ei.get_header_size())

        if len(in_bytes) != self.ei.get_header_size():
            raise ValueError(
                "Header data wrong size: expected %d bytes, got %d" %
                (self.ei.get_header_size(), len(in_bytes))
                )

        hdr = self.ei.process_header(in_bytes)

        img_data = ""
        while len(img_data) < self.ei.get_image_size():
            in_bytes = sock.recv(4096)
            img_data += in_bytes

        if len(img_data) != self.ei.get_image_size():
            if len(img_data) - 8 == self.ei.get_image_size():
                print "seeing byte error, FIXME"
                img_data = img_data[:-8]
            else:
                raise ValueError("Image data wrong size: expected %d bytes, got %d" %
                    (self.ei.get_image_size(), len(img_data)))

        new_ei = self.ei.process_image(img_data)
        if new_ei:
            if isinstance(new_ei, nb.Nifti1Image):
                filename = self.save_nifti(new_ei)
                self.mutex.acquire()
                self.filename_stack.append(filename)
                self.mutex.release()
            else:
                print "ERROR: unknown image type:\n"
                print new_ei
        else:
            self.stop()

    def save_nifti(self, img):
        """Save a received image to a file.
        """

        filename = os.path.join(self.save_location,
                                'img-%05d.nii.gz' % self.save_volume_index)
        img.to_filename(filename)
        self.save_volume_index += 1
        return filename

def parse_args(args):
    """Parse command line arguments.

    USED IN STANDALONE MODE ONLY
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-H", "--host", default="localhost",
                        help="Name of the host to run the image receiver on.")
    parser.add_argument("-p", "--port", type=int, default="15000",
                        help="Port to run the image receiver on.")
    parser.add_argument("-d", "--save_directory", default=".",
                        help="Directory to save images to.")
    return parser.parse_args()

def main(argv):
    """Main entry. Just starts the server and waits for it to finish.

    USED IN STANDALONE MODE ONLY
    """

    args = parse_args(argv)
    receiver = ImageReceiver(args)
    receiver.start()
    while(receiver._is_running):
        sleep(1)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

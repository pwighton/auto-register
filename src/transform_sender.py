"""class to transmit an affine transform over a socket.
"""

import re
import struct

from tcpip_server import ThreadedTCPServer

class TransformSender:

    def __init__(self, port):
        self._port = port

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

    def start(self):
        pass

    def send(self, transform_file):

        if not self.ready:
            print "Transform sender not ready, can't send."
            return False

        total_sent = 0

        transform = TransformSender.read_transform_file(transform_file)
        if len(transform) == 0:
            print "No transform to send"
            return False

        to_send = len(transform)

        while total_sent < to_send:
            sent = self._socket.send(transform[total_sent:])
            if sent == 0:
                print "Socket connection broken"
                break

            total_sent = total_sent + sent

            self._socket.send('\x00')

        return total_sent == to_send

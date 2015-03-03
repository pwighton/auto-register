"""class to transmit an affine transform over a socket.
"""

import re
import struct

from tcpip_server import ThreadedTCPServer

class TransformSender(object):

    def __init__(self, port):
        self._port = port
        self._server = None
        self._transforms_to_send = []

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
        if self._server is not None and self._server.is_running():
            raise RuntimeError('Server already running')

        server = ThreadedTCPServer(("192.168.2.5", self._port), self.process_data)
        ip, port = server.server_address
        print "Transform sender running at %s on port %d" % (ip, port)
        self._server = server

    def send(self, transform_file):

        if not self._server.is_running():
            print "Transform sender not running, can't send."
            return False

        total_sent = 0

        transform = TransformSender.read_transform_file(transform_file)
        if len(transform) == 0:
            print "No transform to send"
            return False

        self._transforms_to_send.append(transform)

        return True

    def process_data(self, sock):
        in_bytes = sock.recv(4096)

        def send_string(data):
            chars_to_send = len(data)

            sent = 0
            while sent < chars_to_send:
                sent_now = sock.send(data[sent:])
                if sent_now <= 0:
                    return False

                sent += sent_now
            return True

        if in_bytes != "ping":
            send_string("unrecognized request %s" % in_bytes)
            return

        if len(self._transforms_to_send) == 0:
            send_string("none\0")
            return

        transform = self._transforms_to_send.pop()
        if not send_string(transform + "\0"):
            print "Error sending transform"
        else:
            print "Transmitted transform string"

    def stop(self):
        if self._server is not None:
            self._server.shutdown()
            self._server = None

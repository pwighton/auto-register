"""class to transmit an affine transform over a socket.
"""

import re
import struct
import time

from tcpip_server import ThreadedTCPServer

class TransformSender(object):

    def __init__(self, host, port, apply_mode):
        self._host = host
        self._port = port
        self._server = None
        self._transforms_to_send = []
        self._apply_mode = apply_mode
        self._delay = 5
        self._delay_timer = time.time()

    def start(self):
        if self._server is not None and self._server.is_running():
            raise RuntimeError('Server already running')

        server = ThreadedTCPServer((self._host, self._port), self.process_data)
        ip, port = server.server_address
        print "Transform sender running at %s on port %d" % (ip, port)
        self._server = server

    def send(self, transform):

        if not self._server.is_running():
            print "Transform sender not running, can't send."
            return False

        total_sent = 0

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

        if ((self._apply_mode and
             self._delay_timer - time.time() < self._delay) or
            len(self._transforms_to_send) == 0):
            send_string("none\0")
            return

        if self._apply_mode:
            transform = self._transforms_to_send[0]
            delay_timer = time.time()
        else:
            transform = self._transforms_to_send.pop()

        if not send_string(transform + "\0"):
            print "Error sending transform"
        else:
            print "Transmitted transform string"



    def stop(self):
        if self._server is not None:
            self._server.shutdown()
            self._server = None

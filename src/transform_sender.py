import re
import struct
import time

from tcpip_server import ThreadedTCPServer

class TransformSender(object):
    """Sends a transformation to an external computer. At the moment, this
    transformation is only suitable for use in the Siemens AutoAlign
    system (or other systems that use a compatible coordinate system).

    This "sender" is actually a TCP server that listens for "ping"
    requests from an external receiver. When a transform is available
    to send, the server responds to the ping request by transmitting
    the transform.

    """

    def __init__(self, host, port):
        """Initialize the networking parameters and internal state.
        """
        self._host = host
        self._port = port
        self._server = None
        self._transforms_to_send = []
        self._state = ""

    def start(self):
        if self._server is not None and self._server.is_running():
            raise RuntimeError('Server already running')

        server = ThreadedTCPServer((self._host, self._port), self.process_data)
        ip, port = server.server_address
        print "Transform sender running at %s on port %d" % (ip, port)
        self._server = server

    def set_state(self, state):
        """Set the internal state of the sender. If the state is anything
        other than an empty string, no transforms will be sent. This
        is useful for disabling sending old transforms while an image
        registration or network transfer is in progress.

        """
        self._state = state

    def clear_state(self):
        """Set the interal state to an empty string (allowing transformations
        to be sent).

        """
        self._state = ""

    def send(self, transform):
        """Queue a transform for sending to an external requestor.
        """

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
        """Callback when a "ping" request is received. This checks that the
        request bytes are actually "ping", and replies with the string
        "none" if there are no transforms queued for sending, and
        otherwise sends the first transform available.

        """

        in_bytes = sock.recv(4096)

        def send_string(data):
            data += "\0"
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
            if self._state == "":
                send_string("none")
            else:
                send_string(self._state)
            return

        transform = self._transforms_to_send.pop()

        if not send_string(transform):
            print "Error sending transform"
        else:
            print "Transmitted transform string"

    def stop(self):
        if self._server is not None:
            self._server.shutdown()
            self._server = None

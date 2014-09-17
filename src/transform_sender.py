"""class to transmit an affine transform over a socket.
"""

import socket
import struct

class TransformSender:

    def __init__(self, host, port):
        self._socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self._socket.connect(host, port)

    def pack754_32(f):
        return pack754(f, 32, 8)

    def pack754_64(f):
        return pack754(f, 64, 11)

    def pack754(f, bits, expbits)
        """Python translation of the code here:
        http://beej.us/guide/bgnet/output/html/singlepage/bgnet.html#serialization
        """
        significandbits = bits - expbits - 1

        if (f == 0.0) return 0

        if f < 0:
            sign = 1
            fnorm = -f
        else:
            sign = 0
            fnorm = f

        shift = 0
        while fnorm >= 2.0:
            fnorm /= 2.0
            shift += 1

        while fnorm < 1.0:
            fnorm *= 2.0
            shift -= 1

        fnorm = fnorm - 1.0

        significand = fnorm * ((1L<<significandbits) + 0.5f)

        exp = shift + ((1<<(expbits-1)) - 1); // shift + bias

        return (sign<<(bits-1)) | (exp<<(bits-expbits-1)) | significand

    def read_transform_file(filename):
        # TODO
        pass

    def pack_transform(transform):
        transform_bytes = []

        for row in xrange(3):
            for col in xrange(4):
                transform_bytes.append(struct.pack('!d', transform[row][col])

        return transform_bytes

    def send(self, transform):

        total_sent = 0

        transform_bytes = pack_transform(transform)
        to_send = len(transform_bytes)

        while total_sent < to_send:
            sent = self._socket.send(transform_bytes[total_sent:])
            if sent == 0:
                print "Socket connection broken"
                break

            total_sent = total_sent + sent

        return total_sent == to_send

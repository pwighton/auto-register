"""Threaded TCP/IP server that will notify a callback on a new connection.
"""

import socket
import SocketServer
import threading

SocketServer.TCPServer.allow_reuse_address = True

class ThreadedTCPRequestHandler(SocketServer.BaseRequestHandler):
    def __init__(self, callback, *args, **keys):
        self.callback = callback
        SocketServer.BaseRequestHandler.__init__(self, *args, **keys)

    def handle(self):
        self.callback(self.request)

def handler_factory(callback):
    def createHandler(*args, **keys):
        return ThreadedTCPRequestHandler(callback,  *args, **keys)
    return createHandler

def process_data_callback(callback, sock):
    callback(sock)

class ThreadedTCPServer(SocketServer.ThreadingMixIn, SocketServer.TCPServer):

    def __init__(self, address, callback):
        SocketServer.TCPServer.__init__(self,
            address, handler_factory(callback))

        # Start a thread with the server -- that thread will then start one
        # more thread for each request
        server_thread = threading.Thread(target=self.serve_forever)
        # Exit the server thread when the main thread terminates
        server_thread.daemon = True
        server_thread.start()
        self._is_running = True

    def is_running(self):
        return self._is_running

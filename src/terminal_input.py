"""class to interact with the user through a terminal.
"""

import fcntl
import os
import sys
import termios
import threading

class TerminalInput(object):

    def __init__(self, disabled):

        self._chars = []
        self._thread = None
        self._mutex = threading.Lock()
        self._should_stop = False

        self._disabled = disabled
        if self._disabled:
            return

        # setup the terminal for interative character capture
        self._fd = sys.stdin.fileno()

        self._oldterm = termios.tcgetattr(self._fd)
        newattr = termios.tcgetattr(self._fd)
        newattr[3] = newattr[3] & ~termios.ICANON & ~termios.ECHO
        termios.tcsetattr(self._fd, termios.TCSANOW, newattr)

        self.oldflags = fcntl.fcntl(self._fd, fcntl.F_GETFL)
        fcntl.fcntl(self._fd, fcntl.F_SETFL, self.oldflags | os.O_NONBLOCK)

    def start(self):
        if self._disabled:
            return

        self._thread = threading.Thread(target=self.run)
        self._thread.start()

    def run(self):
        while not self._should_stop:
            try:
                c = sys.stdin.read(1)

                self._mutex.acquire()
                self._chars.append(c)
                self._mutex.release()
            except IOError:
                pass

    def stop(self):
        if self._disabled:
            return

        self._should_stop = True
        self._thread.join()
        print "Resetting terminal"
        termios.tcsetattr(self._fd, termios.TCSAFLUSH, self._oldterm)
        fcntl.fcntl(self._fd, fcntl.F_SETFL, self.oldflags)

    def get_char(self):

        self._mutex.acquire()
        if len(self._chars) == 0:
            c = None
        else:
            c = self._chars.pop()
        self._mutex.release()

        return c

#ifndef VXWORKS

#include "AffineTransformReceiver.h"

#include "MrCommon\MrNFramework\MrTrace\MPCUTrace\MPCUTrace.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#pragma comment(lib, "ws2_32.lib")

namespace { // variables and functions for internal use

// A TCP/IP server that listens on a socket for data representing an
// affine transform. This function is run in a thread spawned by the
// AffineTransformReceiver::startListening() method.
DWORD WINAPI ListenForTransform(LPVOID lpParam) {

  // TODO: better error reporting

  AffineTransformReceiver *receiver =
    static_cast<AffineTransformReceiver*>(lpParam);

  WSADATA wsaData;
  if (WSAStartup(0x202, &wsaData) == SOCKET_ERROR) {
    WSACleanup();
    return 1;
  }

  sockaddr_in local;
  local.sin_family = AF_INET;
  local.sin_addr.s_addr = INADDR_ANY;
  local.sin_port = htons(receiver->getPortNumber());

  SOCKET listen_socket = socket(AF_INET, SOCK_STREAM, 0);
  if (listen_socket == INVALID_SOCKET){
    WSACleanup();
    return 1;
  }

  if (bind(listen_socket, (sockaddr*) &local, sizeof(local)) == SOCKET_ERROR) {
    WSACleanup();
    return 1;
  }

  if (listen(listen_socket, 5) == SOCKET_ERROR) {
    WSACleanup();
    return 1;
  }

  TRACE_PUT0(TC_ALWAYS, TF_SEQ, "Started listening for transforms");

  const int MAX_TRANSFORM_BYTES = 1024;
  char buffer[MAX_TRANSFORM_BYTES];
  SOCKET msg_sock = 0;

  // receive until signaled by the receiver
  while (receiver->getRunning()) {

    if (!msg_sock) {
      msg_sock = accept(listen_socket, NULL, NULL);
      if (msg_sock == INVALID_SOCKET) {
        WSACleanup();
        return 1;
      }
    }

    int num_bytes = 0;
    bool success = false;
    for (; num_bytes < MAX_TRANSFORM_BYTES; num_bytes++) {
      int received = recv(msg_sock, buffer + num_bytes, 1, 0);

      if (received != 1) {
        TRACE_PUT1(TC_ALWAYS, TF_SEQ, "Error receiving on the socket, code: %d", received);
        closesocket(msg_sock);
        msg_sock = 0;
        break;
      }

      if (*(buffer + num_bytes) == '\0') {
        success = true;
        break;
      }
    }

    if (success) {
      receiver->setTransformFromString(buffer);
      // send the acknowledgement
      const char *ack = "received";
      send(msg_sock, ack, strlen(ack), 0);
    }
    else {
      TRACE_PUT0(TC_ALWAYS, TF_SEQ, "Unknown error receiving transform");
      closesocket(msg_sock);
      msg_sock = 0;
    }
  }

  if (msg_sock != 0) {
    closesocket(msg_sock);
  }

  closesocket(listen_socket);

  receiver->setThreadExit(true);

  return 0;
}


} // anonymouse namespace

AffineTransformReceiver::AffineTransformReceiver(int port_number)
  : port_number(port_number)
  , should_stop(false)
  , thread_handle(0)
  , thread_exit(false)
  , thread_exit_mutex(0)
  , has_transform(false)
  , has_transform_mutex(0)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      matrix[i][j] = 0.0;
    }
  }
}

AffineTransformReceiver::~AffineTransformReceiver() {
  if (thread_handle != 0) {
    stopListening();
  }
}

bool AffineTransformReceiver::startListening() {
  if (thread_handle != 0) {
    return false;
  }

  thread_exit_mutex = CreateMutex(NULL, false, NULL);
  thread_exit = false;

  thread_handle = CreateThread(NULL, 0, ListenForTransform,
    static_cast<void*>(this), 0, NULL);

  return true;
}

bool AffineTransformReceiver::stopListening() {

  should_stop = true;

  while (!getThreadExit()) {
    Sleep(10);
  }

  CloseHandle(thread_handle);
  thread_handle = 0;

  CloseHandle(thread_exit_mutex);
  thread_exit_mutex = 0;

  return true;
}

bool AffineTransformReceiver::getThreadExit() const {
  bool ret = false;
  WaitForSingleObject(thread_exit_mutex, INFINITE);
  ret = thread_exit;
  ReleaseMutex(thread_exit_mutex);
  return ret;
}

void AffineTransformReceiver::setThreadExit(bool value) {
  WaitForSingleObject(thread_exit_mutex, INFINITE);
  thread_exit = value;
  ReleaseMutex(thread_exit_mutex);
}

bool AffineTransformReceiver::getHasTransform() const {
  bool ret = false;
  WaitForSingleObject(has_transform_mutex, INFINITE);
  ret = has_transform;
  ReleaseMutex(has_transform_mutex);
  return ret;
}

void AffineTransformReceiver::setHasTransform(bool value) {
  WaitForSingleObject(has_transform_mutex, INFINITE);
  has_transform = value;
  ReleaseMutex(has_transform_mutex);
}

double AffineTransformReceiver::getTransformMatrixEl(int r, int c) const {
  double ret = 0.0;
  WaitForSingleObject(has_transform_mutex, INFINITE);
  ret = matrix[r][c];
  ReleaseMutex(has_transform_mutex);
  return ret;
}

bool AffineTransformReceiver::setTransformFromString(char *string) {

  TRACE_PUT1(TC_ALWAYS, TF_SEQ, "%s", string);

  WaitForSingleObject(has_transform_mutex, INFINITE);
  int start = 0;
  int end = 0;
  int row = 0;
  int col = 0;
  while (1) {
    while (!isspace(string[end])) {
      end++;
    }

    char *double_chars = new char[end - start + 1];
    strncpy(double_chars, string + start, end - start);
    double_chars[end - start] = '\0';
    matrix[row][col] = atof(double_chars);

    if (++col >= 4) {
      col = 0;
      if (++row >= 3) {
        break;
      }
    }

    start = end + 1;
    end = start;
  }
  ReleaseMutex(has_transform_mutex);
  setHasTransform(true);
  stopListening();

  return true;
}

#endif

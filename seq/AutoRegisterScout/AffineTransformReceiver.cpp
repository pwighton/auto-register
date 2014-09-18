#ifndef VXWORKS

#include "AffineTransformReceiver.h"

#include "MrCommon\MrNFramework\MrTrace\MPCUTrace\MPCUTrace.h"

#pragma comment(lib, "ws2_32.lib")

using std::string;

AffineTransformReceiver::AffineTransformReceiver(const string &host,
                                                 int port_number)
  : host(host)
  , port_number(port_number)
  , ready_for_network(false)
{
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      matrix[i][j] = 0.0;
    }
  }
}

AffineTransformReceiver::~AffineTransformReceiver() {
  if (ready_for_network) {
    WSACleanup();
    ready_for_network = false;
  }
}

bool AffineTransformReceiver::initNetwork() {
  WSADATA wsaData;
  if (WSAStartup(0x202, &wsaData) == SOCKET_ERROR) {
    WSACleanup();
    ready_for_network = false;
    return false;
  }

  ready_for_network = true;
  return true;
}

bool AffineTransformReceiver::checkForTransform() {

  if (!ready_for_network) {
    if (!initNetwork()) {
      return false;
    }
  }

  sockaddr_in remote;
  remote.sin_family = AF_INET;
  remote.sin_addr.s_addr = inet_addr(host.c_str());
  remote.sin_port = htons(port_number);

  SOCKET transform_socket = socket(AF_INET, SOCK_STREAM, 0);
  if (transform_socket == INVALID_SOCKET){
    return false;
  }

  TRACE_PUT0(TC_ALWAYS, TF_SEQ, "Checking for transform");

  if (connect(transform_socket, (SOCKADDR*) &remote, sizeof(remote) ==
              SOCKET_ERROR)) {
    return false;
  }

  // send the request for a transform
  const string request("ping");
  int sent = send(transform_socket, request.c_str(), request.length(), 0);
  if (sent != request.length()) {
    TRACE_PUT0(TC_ALWAYS, TF_SEQ, "Failed to send request.");
    return false;
  }

  // receive the response
  const int MAX_TRANSFORM_BYTES = 1024;
  const string NO_TRANSFORM_MESSAGE("none");
  char buffer[MAX_TRANSFORM_BYTES];

  int num_bytes = 0;
  bool success = false;
  for (; num_bytes < MAX_TRANSFORM_BYTES; num_bytes++) {
    int received = recv(msg_sock, buffer + num_bytes, 1, 0);

    if (received != 1) {
      TRACE_PUT1(TC_ALWAYS, TF_SEQ, "Error receiving on the socket, code: %d",
                 received);
      return false;
    }

    if (*(buffer + num_bytes) == '\0') {
      success = true;
      break;
    }
  }

  if (!success && num_bytes == MAX_TRANSFORM_BYTES) {
    TRACE_PUT0(TC_ALWAYS, TF_SEQ,
               "Error: transform data too long or data corruption.");
    return false;
  }
  else if (!success) {
    TRACE_PUT0(TC_ALWAYS, TF_SEQ, "Unknown error receiving transform.");
    return false;
  }
  else if (NO_TRANSFORM_MESSAGE == buffer) {
    return false;
  }

  setTransformFromString(buffer);
  return true;
}

bool AffineTransformReceiver::setTransformFromString(const string &str) {

  TRACE_PUT1(TC_ALWAYS, TF_SEQ, "%s", string);

  int start = 0;
  int end = 0;
  int row = 0;
  int col = 0;
  while (1) {
    while (!isspace(str[end])) {
      end++;
    }

    char *double_chars = new char[end - start + 1];
    strncpy(double_chars, str + start, end - start);
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

  return true;
}

#endif

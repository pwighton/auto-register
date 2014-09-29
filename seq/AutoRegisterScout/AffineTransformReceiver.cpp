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

  hostent *host_info = NULL;
  if (isalpha(host[0])) {
    host_info = gethostbyname(host.c_str());
  }
  else  {
    unsigned int addr = inet_addr(host.c_str());
    host_info = gethostbyaddr((char *) &addr, 4, AF_INET);
  }

  if (host_info == NULL) {
    TRACE_PUT0(TC_ALWAYS, TF_SEQ, "Error resolving host");
  }

  sockaddr_in remote;
  memset(&remote, 0, sizeof(remote));
  memcpy(&(remote.sin_addr), host_info->h_addr, host_info->h_length);
  remote.sin_family = host_info->h_addrtype;
  remote.sin_port = htons(port_number);

  SOCKET transform_socket = socket(AF_INET, SOCK_STREAM, 0);
  if (transform_socket == INVALID_SOCKET){
    TRACE_PUT0(TC_ALWAYS, TF_SEQ, "Error allocating socket");
    return false;
  }

  TRACE_PUT0(TC_ALWAYS, TF_SEQ, "Checking for transform");

  if (connect(transform_socket, (SOCKADDR*) &remote, sizeof(remote)) ==
              SOCKET_ERROR) {
    TRACE_PUT0(TC_ALWAYS, TF_SEQ, "Failed to connect to server.");
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

  int received = recv(transform_socket, buffer, MAX_TRANSFORM_BYTES, 0);

  if (received <= 0) {
    TRACE_PUT1(TC_ALWAYS, TF_SEQ,
               "Error receiving on the socket, code: ", received);
    return false;
  }

  TRACE_PUT1(TC_ALWAYS, TF_SEQ, "Response: %s.", buffer);

  if (received == MAX_TRANSFORM_BYTES) {
    TRACE_PUT0(TC_ALWAYS, TF_SEQ,
               "Error: transform data too long or data corruption.");
    return false;
  }
  else if (NO_TRANSFORM_MESSAGE == buffer) {
    return false;
  }

  setTransformFromString(buffer);
  return true;
}

bool AffineTransformReceiver::setTransformFromString(const char *str) {

  TRACE_PUT1(TC_ALWAYS, TF_SEQ, "%s", str);

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

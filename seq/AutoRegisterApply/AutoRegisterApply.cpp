#ifndef VXWORKS

#include "AutoRegisterApply.h"

#include "MrServers/MrProtSrv/MrProt/MrProt.h"
#include "MrServers/MrProtSrv/MrProt/MrVector.h"

#include "MrCommon/MrNFramework/MrTrace/MPCUTrace/MPCUTrace.h"

#include "MrServers/MrImaging/seq/AutoRegisterApply/mgh_isometry/MGH_AffineTransform.h"

#define WIN32_LEAN_AND_MEAN
#include<winsock2.h>
#pragma comment(lib, "ws2_32.lib")

using std::string;

namespace {

const string DEFAULT_HOST = "192.168.1.2";
const int DEFAULT_PORT = 15001;

// Local class to support retreiving transform data from an external
// computer. Each time checkForTransform() is called, a TCP/IP
// connection is made to the specified host to check for a new
// transform.
class AffineTransformReceiver {
 public:
  AffineTransformReceiver(const std::string &host, int port_number)
    : host(host)
    , port_number(port_number)
    {}

  ~AffineTransformReceiver() {}

  // check whether a new transform is available.
  bool checkForTransform() {
    WSADATA wsaData;
    if (WSAStartup(0x202, &wsaData) == SOCKET_ERROR) {
      WSACleanup();
      return false;
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
      return false;
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

    WSACleanup();

    return true;
  }

  // retrieve a single affine matrix element at [r,c].
  // no bounds checking or validation is performed.
  double getTransformMatrixEl(int r, int c) const {
    TRACE_PUT3(TC_ALWAYS, TF_SEQ, "%d %d %f", r, c, transform.matrix[r][c]);
    return transform.matrix[r][c];
  }

  const LIB_NAMESPACE::AffineTransform& getTransformMatrix() const {
    return transform;
  }

 private:

  // setup the transform matrix based on a string of doubles.
  bool setTransformFromString(const char *str) {

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
      transform.matrix[row][col] = atof(double_chars);

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

  std::string host;
  int port_number;

  LIB_NAMESPACE::AffineTransform transform;
};

} // anonymous namespace

AutoRegisterApply::AutoRegisterApply()
  : host(DEFAULT_HOST)
  , port_number(DEFAULT_PORT)
{}

AutoRegisterApply::AutoRegisterApply(const string &host, int port_number)
  : host(host)
  , port_number(port_number)
{}

AutoRegisterApply::~AutoRegisterApply()
{}

bool AutoRegisterApply::applyToProtocol(MrProt *pMrProt) const {
  // check for transform
  AffineTransformReceiver receiver(host, port_number);

  if (!receiver.checkForTransform()) {
    TRACE_PUT0(TC_ALWAYS, TF_SEQ, "No transform received");
    return false;
  }

  // apply the transform to the protocol

  // set the axis angle
  double angle, nx, ny, nz;
  receiver.getTransformMatrix().toAxisAngle(&angle, &nx, &ny, &nz);
  TRACE_PUT4(TC_ALWAYS, TF_SEQ, "axis angle: %g %g %g %g", angle, nx, ny, nz);

  pMrProt->sliceGroupList()[0].rotationAngle(angle);
  pMrProt->sliceGroupList()[0].normal(nx, ny, nz);

  // // set the position
  VectorPat<double> pos;
  pos.dSag = receiver.getTransformMatrixEl(0, 3);
  pos.dCor = receiver.getTransformMatrixEl(1, 3);
  pos.dTra = receiver.getTransformMatrixEl(2, 3);
  pMrProt->sliceGroupList()[0].position(pos);

  return true;
}

#endif

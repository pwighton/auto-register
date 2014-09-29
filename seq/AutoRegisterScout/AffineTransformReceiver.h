// Functions to support listening for transform data from an external
// computer. Each time checkForTransform() is called, a TCP/IP
// connection is made to the specified host to check for a new
// transform.

#ifndef VXWORKS

#pragma once

#include <string>

#define WIN32_LEAN_AND_MEAN
#include<winsock2.h>

class AffineTransformReceiver {
 public:
  AffineTransformReceiver(const std::string &host, int port_number);

  ~AffineTransformReceiver();

  // initialize network operations
  bool initNetwork();

  // check whether a new transform is available.
  bool checkForTransform();

  // retrieve a single affine matrix element at [r,c].
  // no bounds checking or validation is performed.
  double getTransformMatrixEl(int r, int c) const {
    return matrix[r][c];
  }

 private:

  // setup the transform matrix based on a string of doubles.
  bool setTransformFromString(const char *);

  bool ready_for_network;

  std::string host;
  int port_number;

  double matrix[3][4];
};

#endif

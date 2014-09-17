// Functions to support listening for transform data from an external
// computer. A TCP/IP server is run in a thread, and sets flags when a
// transform is received. 

#ifndef VXWORKS

#pragma once

#define WIN32_LEAN_AND_MEAN
#include<winsock2.h>

class AffineTransformReceiver {
 public:
  AffineTransformReceiver(int port_number);

  ~AffineTransformReceiver();

  int getPortNumber() const {
    return port_number;
  }

  // initialize and start the server to listen for incoming transforms
  bool startListening();

  // tear down the server and join the listener thread. 
  // is called internally after a transform is received. 
  // should also be called before running the sequence.
  bool stopListening();

  // whether we've been told to stop running
  bool getRunning() const {
    return !should_stop;
  }

  // whether the listener thread has terminated
  // Note: blocks on thread_exit_mutex
  bool getThreadExit() const;

  // the listener thread let's us know it's done by calling this
  // function with 'true'
  // Note: blocks on thread_exit_mutex
  void setThreadExit(bool exited);

  // whether we've received a transform via our listener thread
  // Note: blocks on has_transform_mutex
  bool getHasTransform() const;

  // the listener thread let's us know it's received a transform by
  // calling this with 'true'
  // Note: blocks on has_transform_mutex
  void setHasTransform(bool value);

  // retrieve a single affine matrix element at [r,c]
  // Note: blocks on has_transform_mutex
  double getTransformMatrixEl(int r, int c) const;

  // setup the transform matrix based on a string of doubles.
  // Note: blocks on has_transform_mutex
  bool setTransformFromString(char *string);

 private:

  int port_number;

  double matrix[3][4];
  bool should_stop;

  HANDLE thread_handle;
  bool thread_exit;
  HANDLE thread_exit_mutex;

  bool has_transform;
  HANDLE has_transform_mutex;
};

#endif

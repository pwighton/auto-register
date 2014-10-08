// Methods for applying a transformation found using the AutoRegister
// scout to the slice prescription of another sequence.
//
// USAGE:
//
// Sequence code:
//
// The first time pSeqLim->isContextNormal() is true in
// fSEQPrep(), instantiate an AutoRegisterApply object, then call
// applyToProtocol() on it. This will check with the TCPIP server
// running on the host and port supplied, and if a transformation is
// found, will modify the slice presciption of the passed protocol
// accordingly.
//
// Build environment:
//
// Add to <sequence>.dsp:
//
// In the sources section:
// # Begin Source File
// SOURCE=..\..\AutoRegisterApply\AutoRegisterApply.cpp
// # End Source File
// # Begin Source File
// SOURCE=..\..\AutoRegisterApply\mgh_isometry\MGH_AffineTransform.cpp
// # End Source File
// # Begin Source File
// SOURCE=..\..\AutoRegisterApply\mgh_isometry\MGH_Quaternion.cpp
// # End Source File
//
// In the headers section:
//
// # Begin Source File
// SOURCE=..\..\AutoRegisterApply\AutoRegisterApply.h
// # End Source File
// # Begin Source File
// SOURCE=..\..\AutoRegisterApply\mgh_isometry\MGH_AffineTransform.h
// # End Source File
// # Begin Source File
// SOURCE=..\..\AutoRegisterApply\mgh_isometry\MGH_Quaternion.h
// # End Source File
//
// Add to <sequence>.mk:
// SOURCES += ${IDEA_BASE}/n4/pkg/MrServers/MrImaging/seq/AutoRegisterApply/AutoRegisterApply.cpp
// SOURCES += ${IDEA_BASE}/n4/pkg/MrServers/MrImaging/seq/AutoRegisterApply/mgh_isometry/MGH_AffineTransform.cpp
// SOURCES += ${IDEA_BASE}/n4/pkg/MrServers/MrImaging/seq/AutoRegisterApply/mgh_isometry/MGH_Quaternion.cpp
//
// Add to makefile.trs:
// CPPSOURCES (..\AutoRegisterApply\AutoRegisterApply)
// CPPSOURCES (..\AutoRegisterApply\mgh_isometry\MGH_AffineTransform)
// CPPSOURCES (..\AutoRegisterApply\mgh_isometry\MGH_Quaternion)


#ifndef VXWORKS
#pragma once

#include <string>

class MrProt;

class AutoRegisterApply {

 public:

  AutoRegisterApply();

  AutoRegisterApply(const std::string &host, int port);

  ~AutoRegisterApply();

  // polls a TCPIP server at host:port for an affine transform. If
  // found, the transform is applied to the slice prescription in
  // pMrProt. Returns true if a transform was applied.
  bool applyToProtocol(MrProt *pMrProt) const;

 private:

  std::string host;
  int port_number;

};

#endif

// Utility to serve a nifti image file simulating what would be sent
// by VSend from the scanner.
//
// Modified from servenii4d.cpp from murfi2.

#include<iostream>

#include"ace/SOCK_Stream.h"
#include"ace/SOCK_Connector.h"
#include"ace/SOCK_Stream.h"

#include"nifti1_io.h"

#include"RtExternalSenderImageInfo.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;

void usage(char *execname) {
}

const int DEFAULT_PORT = 15000;

class VSendNifti {
 public:
  VSendNifti(const string &host)
    : address(DEFAULT_PORT, host.c_str())
    , image(NULL) {}

  bool init(const string &nifti_filename) {

    // load the image
    image = nifti_image_read(nifti_filename.c_str(), 0);
    if(image == NULL) {
      cerr << "could not open " << nifti_filename
           << " for reading a nifti image" << endl;
      return false;
    }

    // validate
    if(image->datatype != DT_SIGNED_SHORT) {
      cerr << "only signed short nifti files are supported" << endl;
      return false;
    }

    // build the image info
    strcpy(image_info.magic, "SIMU");
    strcpy(image_info.imageType, "3D");
    strcpy(image_info.scanType, "FLASH");
    strcpy(image_info.dataType, "int16_t");

    image_info.isLittleEndian = true;
    image_info.isMosaic = false;

    image_info.numPixelsRead = image->dim[1];
    image_info.numPixelsPhase = image->dim[2];
    image_info.numSlices = image->dim[3];
    image_info.totalTR = image->dim[4];

    image_info.pixelSpacingReadMM = image->pixdim[1];
    image_info.pixelSpacingPhaseMM = image->pixdim[2];
    image_info.pixelSpacingSliceMM = image->pixdim[3];
    image_info.repetitionTimeMS = 1000 * image->pixdim[4];

    for (int r = 0; r < 3; r++) {
      for (int c = 0; c < 3; c++) {
        image_info.voxelToWorldMatrix[r][c] = image->qto_xyz.m[r][c];
      }
    }

    image_info.voxelToWorldMatrix[0][3] = image->qto_xyz.m[0][3];
    image_info.voxelToWorldMatrix[1][3] = image->qto_xyz.m[1][3];
    image_info.voxelToWorldMatrix[2][3] = image->qto_xyz.m[2][3];

    image_info.isMotionCorrected = false;

    // read the image data
    nifti_image_load(image);

    cout << "read header and data from " << nifti_filename << endl;

    return true;
  }

  ~VSendNifti() {
    free(image);
  };

  void run() {

    size_t total_images = image->dim[4];
    size_t cur_image = 0;

    // Keep making new connections until we've sent the whole series.
    for(; cur_image < total_images &&
          !connector.connect(stream, address); cur_image++) {
      cout << "made connection, sending image " << cur_image + 1 << endl;

      image_info.currentTR = cur_image + 1;

      // send the image header
      stream.send_n(reinterpret_cast<char*>(&image_info),
                    image_info.getHeaderSize());

      // send the image
      short *img_data = (short*) image->data +
        cur_image * image->dim[1] * image->dim[2] * image->dim[3];

      int sent = stream.send_n(img_data, image_info.getDataSize());

      cout << "sent " << sent << endl;

      stream.close();

      if (cur_image < total_images - 1) {
        usleep(1000000 * image->pixdim[4]);
      }
    }

    if (cur_image == 0 && total_images > 0) {
      cerr << "Failed to connect" << endl;
    }
  }

 private:

  // ACE objects for socket connections
  ACE_INET_Addr address;
  ACE_SOCK_Stream stream;
  ACE_SOCK_Connector connector;

  // image
  nifti_image *image;

  // image header to send
  RtExternalImageInfo image_info;
};

int ACE_TMAIN (int argc, ACE_TCHAR *argv[]) {
  if (argc < 2) {
    cerr << "filename is required" << endl;
    cerr << "usage: " << argv[0] << " niftifile " << endl;
    return 1;
  }

  string host = "";
  if (argc >= 3) {
    host = argv[2];
  }

  VSendNifti sender(host);
  if (!sender.init(argv[1])) {
    cerr << "usage: " << argv[0] << " niftifile [host]" << endl;
    return 1;
  }

  sender.run();

  return 0;
}

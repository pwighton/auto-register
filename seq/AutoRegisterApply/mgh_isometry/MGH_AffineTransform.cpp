#include "MGH_AffineTransform.h"

#include <cmath>

#include "MGH_Quaternion.h"

namespace LIB_NAMESPACE {
const double AffineTransform::identityRotation[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
const double AffineTransform::identityTranslation[3] = {0,0,0};


AffineTransform::AffineTransform() {
}

AffineTransform::AffineTransform(const double rotation[3][3], const double translation[3]) {
	int i, j;

	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			matrix[i][j] = rotation[i][j];
		}
	}

	for(i = 0; i < 3; i++) {
		matrix[i][3] = translation[i];
	}

	for(j = 0; j < 3; j++) {
		matrix[3][j] = 0;
	}

	matrix[3][3] = 1;
}


AffineTransform& AffineTransform::operator=(const AffineTransform &other) {
	int i, j;

	for(i = 0; i < 4; i++) {
		for(j = 0; j < 4; j++) {
			matrix[i][j] = other.matrix[i][j];
		}
	}

	return *this;
}

AffineTransform AffineTransform::operator*(const AffineTransform &other) const {
	AffineTransform out;
	// compute values for first three rows (last row is known a priori if these are indeed affine transforms)
	for(int outRow = 0; outRow < 3; outRow++) {
		for(int outCol = 0; outCol < 4; outCol++) {
			out.matrix[outRow][outCol] = matrix[outRow][0] * other.matrix[0][outCol] + matrix[outRow][1] * other.matrix[1][outCol] + matrix[outRow][2] * other.matrix[2][outCol] + matrix[outRow][3] * other.matrix[3][outCol];
		}
	}

	// zero first three columns of last row
	for(int j = 0; j < 3; j++) {
		out.matrix[3][j] = 0;
	}

	// set last element to one
	out.matrix[3][3] = 1;

	return out;
}

AffineTransform AffineTransform::inverse() const {
	double	invRotation [3][3];

	// Dylan Dec 5/13 -- transpose rotation part
	invRotation[0][0] = matrix[0][0];
	invRotation[1][0] = matrix[0][1];
	invRotation[2][0] = matrix[0][2];

	invRotation[0][1] = matrix[1][0];
	invRotation[1][1] = matrix[1][1];
	invRotation[2][1] = matrix[1][2];

	invRotation[0][2] = matrix[2][0];
	invRotation[1][2] = matrix[2][1];
	invRotation[2][2] = matrix[2][2];

	// Dylan Dec 5/13 -- negate translation part
	double invTranslation[3] = {-matrix[0][3], -matrix[1][3], -matrix[2][3]};

	// Dylan Dec 5/13 -- put them together in the opposite order to generate the inverse transform
	return AffineTransform(invRotation, AffineTransform::identityTranslation) * AffineTransform(AffineTransform::identityRotation, invTranslation);
}

void AffineTransform::toAxisAngle(double* angle, double* x, double* y, double* z) const {
  double rot_matrix[3][3];
  RotationPartMatrix(rot_matrix);

  // Convert using the Quaternion representation's method.
  Quaternion::fromRotationMatrix(rot_matrix).toAxisAngle(angle, x, y, z);
}

std::ostream& operator<<(std::ostream& ost, const AffineTransform& toOut) {
    ost << toOut.matrix[0][0] << " | " << toOut.matrix[0][1] << " | " << toOut.matrix[0][2] << " | " << toOut.matrix[0][3] << std::endl;
    ost << toOut.matrix[1][0] << " | " << toOut.matrix[1][1] << " | " << toOut.matrix[1][2] << " | " << toOut.matrix[1][3] << std::endl;
    ost << toOut.matrix[2][0] << " | " << toOut.matrix[2][1] << " | " << toOut.matrix[2][2] << " | " << toOut.matrix[2][3] << std::endl;
    ost << toOut.matrix[3][0] << " | " << toOut.matrix[3][1] << " | " << toOut.matrix[3][2] << " | " << toOut.matrix[3][3] << std::endl;
	return ost;
}


// Dylan May 6/13 - This is okay away from points where gimbal lock becomes an issue. We shouldn't be working with big rotations in Euler angles anyway, so this is fine.
void AffineTransform::RotationPartZYXEulerAngles(double *angles) const {
	angles[0] = atan2(matrix[1][0], matrix[0][0]);
	angles[1] = -asin(matrix[2][0]);
	angles[2] = atan2(matrix[2][1], matrix[2][2]);
}

void AffineTransform::RotationPartMatrix(double rot_matrix[3][3]) const {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rot_matrix[i][j] = matrix[i][j];
    }
  }
}

}

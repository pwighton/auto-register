#include "MGH_Quaternion.h"

// Dylan Dec 19/12
// For a sensible introduction to dual quaternions (with references), check out:
// Dual Quaternions for Rigid Transformation Blending
// Kavan et al.
// https://www.cs.tcd.ie/publications/tech-reports/reports.06/TCD-CS-2006-46.pdf
// Also, lots of material on representing position with dual quaternions in the
// robotics literature.

#include <cmath>

namespace LIB_NAMESPACE {

const double Quaternion::pi = 3.14159265358979323846;
const double Quaternion::pix2 = 3.14159265358979323846 * 2.0;

Quaternion::Quaternion() : w(1), x(0), y(0), z(0) {
}

Quaternion::Quaternion(double initW, double initX, double initY, double initZ) : w(initW), x(initX), y(initY), z(initZ) {
}

double Quaternion::norm() const {
	return sqrt(w * w + x * x + y * y + z * z);
}


Quaternion Quaternion::normalize() const {
	double normScale = 1.0 / norm();

	return Quaternion(w * normScale, x * normScale, y * normScale, z * normScale);
}

Quaternion Quaternion::conjugate() const{
	return Quaternion(w, -x, -y, -z);
}

void Quaternion::toAxisAngle(double *radians, double *axisX, double *axisY, double *axisZ) {
	*radians = 2.0 * acos(w);
	double invSHalfAngle = 1.0 / sin(*radians / 2.0);
	*axisX = x * invSHalfAngle;
	*axisY = y * invSHalfAngle;
	*axisZ = z * invSHalfAngle;

	while(*radians > Quaternion::pi) *radians -= pix2;
	while(*radians < - Quaternion::pi) *radians += pix2;
}

Quaternion Quaternion::fromAxisAngle(double radians, double ax, double ay, double az) {
  double w = cos(radians / 2);
  double x = ax * sin(radians / 2);
  double y = ay * sin(radians / 2);
  double z = az * sin(radians / 2);
  return Quaternion(w, x, y, z);
}

void Quaternion::toRotationMatrix(double rotMatrix[3][3]) {
	double wSq = w * w;
	double xSq = x * x;
	double ySq = y * y;
	double zSq = z * z;

	double wx2 = w*x*2.0;
	double wy2 = w*y*2.0;
	double wz2 = w*z*2.0;

	double xy2 = x*y*2.0;
	double xz2 = x*z*2.0;

	double yz2 = y*z*2.0;

	rotMatrix[0][0] = wSq + xSq - ySq - zSq;
	rotMatrix[0][1] = xy2 - wz2;
	rotMatrix[0][2] = xz2 + wy2;

	rotMatrix[1][0] = xy2 + wz2;
	rotMatrix[1][1] = wSq - xSq + ySq - zSq;
	rotMatrix[1][2] = yz2 - wx2;

	rotMatrix[2][0] = xz2 - wy2;
	rotMatrix[2][1] = yz2 + wx2;
	rotMatrix[2][2] = wSq - xSq - ySq + zSq;
}

Quaternion Quaternion::fromRotationMatrix(const double matrix[3][3]) {

	// Dylan Dec 19/12
	// This algorithm taken from http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm

	// Also tried the algorithm in
	//		Animating rotation with quaternion curves.
	//		Ken Shoemake
	//		Computer Graphics 19(3):245-254,  1985
	//		http://portal.acm.org/citation.cfm?doid=325334.325242
	// but you'll find that it's not numerically stable without some truncation epsilon.
	// The algorithm we're using now doesn't require us to pick some arbitrary epsilon, so
	// I like it better.

	double tr = matrix[0][0] + matrix[1][1] + matrix[2][2];
	double w, x, y, z, S, SInv;

	if (tr > 0) {
	  S = sqrt(tr+1.0) * 2; // S=4*qw
	  SInv = 1.0 / S;
	  w = 0.25 * S;
	  x = (matrix[2][1] - matrix[1][2]) * SInv;
	  y = (matrix[0][2] - matrix[2][0]) * SInv;
	  z = (matrix[1][0] - matrix[0][1]) * SInv;
	} else if ((matrix[0][0] > matrix[1][1]) && (matrix[0][0] > matrix[2][2])) {
	  S = sqrt(1.0 + matrix[0][0] - matrix[1][1] - matrix[2][2]) * 2; // S=4*qx
	  SInv = 1.0 / S;
	  w = (matrix[2][1] - matrix[1][2]) * SInv;
	  x = 0.25 * S;
	  y = (matrix[0][1] + matrix[1][0]) * SInv;
	  z = (matrix[0][2] + matrix[2][0]) * SInv;
	} else if (matrix[1][1] > matrix[2][2]) {
	  S = sqrt(1.0 + matrix[1][1] - matrix[0][0] - matrix[2][2]) * 2; // S=4*qy
	  SInv = 1.0 / S;
	  w = (matrix[0][2] - matrix[2][0]) * SInv;
	  x = (matrix[0][1] + matrix[1][0]) * SInv;
	  y = 0.25 * S;
	  z = (matrix[1][2] + matrix[2][1]) * SInv;
	} else {
	  S = sqrt(1.0 + matrix[2][2] - matrix[0][0] - matrix[1][1]) * 2; // S=4*qz
	  SInv = 1.0 / S;
	  w = (matrix[1][0] - matrix[0][1]) * SInv;
	  x = (matrix[0][2] + matrix[2][0]) * SInv;
	  y = (matrix[1][2] + matrix[2][1]) * SInv;
	  z = 0.25 * S;
	}

	return Quaternion(w,x,y,z);
}

Quaternion& Quaternion::operator=(const Quaternion &other) {
	w = other.w;
	x = other.x;
	y = other.y;
	z = other.z;
	return *this;
}

Quaternion Quaternion::operator*(const Quaternion &other) const {
	return Quaternion(
		w * other.w - x * other.x - y * other.y - z * other.z,
		w * other.x + x * other.w + y * other.z - z * other.y,
		w * other.y - x * other.z + y * other.w + z * other.x,
		w * other.z + x * other.y - y * other.x + z * other.w
		);
}

Quaternion Quaternion::operator+(const Quaternion &other) const {
	return Quaternion(
		w + other.w,
		x + other.x,
		y + other.y,
		z + other.z
		);
}

Quaternion Quaternion::operator-(const Quaternion &other) const {
	return Quaternion(
		w - other.w,
		x - other.x,
		y - other.y,
		z - other.z
		);
}

std::ostream& operator<<(std::ostream& ost, const Quaternion& toOut) {
	ost << toOut.w << ", " << toOut.x << ", " << toOut.y << ", " << toOut.z;
	return ost;
}



DualQuaternion::DualQuaternion() : real(1,0,0,0), dual(0,0,0,0) {
}


DualQuaternion::DualQuaternion(const Quaternion initReal, const Quaternion initDual) : real(initReal), dual(initDual) {
}

DualQuaternion DualQuaternion::fromRotation(const double radians, const double x, const double y, const double z) {
	double halfAngle = radians * 0.5;
	double sHalfAngle = sin(halfAngle);
	return DualQuaternion(Quaternion(cos(halfAngle), x * sHalfAngle, y * sHalfAngle, z * sHalfAngle), Quaternion(0.0, 0.0, 0.0, 0.0));
}

// Dylan Dec 10/12 -- VC++ doesn't have copysign in all versions, so we'll just make a platform-independent one here
double localCopySign(double x, double y) {
	return fabs(x) * ((y > 0) - (y < 0));
}


DualQuaternion DualQuaternion::fromRotationMatrix(const double matrix[3][3]) {
  return DualQuaternion(Quaternion::fromRotationMatrix(matrix), Quaternion(0,0,0,0));
}

/*
DualQuaternion DualQuaternion::fromRotationMatrix(const double matrix[3][3]) {
	double xx = matrix[0][0];
	double yy = matrix[1][1];
	double zz = matrix[2][2];

	double w, x, y, z;

	// Dylan Dec 10/12 -- This is taken from:
	//		Animating rotation with quaternion curves.
	//		Ken Shoemake
	//		Computer Graphics 19(3):245-254,  1985
	//		http://portal.acm.org/citation.cfm?doid=325334.325242

	double wSq = 0.25 * (xx + yy + zz + 1.0);

	if(wSq > 0) {
		w = sqrt(wSq);
		double invW4 = 1.0 / (4.0 * w);

		x = (matrix[1][2] - matrix[2][1]) * invW4;
		y = (matrix[2][0] - matrix[0][2]) * invW4;
		z = (matrix[0][1] - matrix[1][0]) * invW4;
	}
	else {
		w = 0;
		double xSq = -0.5 * (yy + zz);

		if(xSq > 0) {
			x = sqrt(xSq);
			double invX2 = 1.0 / (2.0 * x);

			y = matrix[0][1] * invX2;
			z = matrix[0][2] * invX2;
		}
		else {
			x = 0;

			double ySq = 0.5 * (1.0 - zz);

			if(ySq > 0) {
				y = sqrt(ySq);
				z = matrix[1][2] / (2.0 * y);
			}
			else {
				y = 0;
				z = 1;
			}
		}
	}

	return DualQuaternion(Quaternion(w,x,y,z), Quaternion(0,0,0,0));
}
*/

DualQuaternion DualQuaternion::fromTranslation(const double x, const double y, const double z) {
	return DualQuaternion(Quaternion(1.0, 0.0, 0.0, 0.0), Quaternion(0.0, x * 0.5, y * 0.5, z * 0.5));
}

DualQuaternion & DualQuaternion::operator=(const DualQuaternion &other) {
	real = other.real;
	dual = other.dual;
	return *this;
}

DualQuaternion DualQuaternion::operator*(const DualQuaternion &other) const {
	return DualQuaternion(
		real * other.real,
		(real * other.dual) + (dual * other.real)
		);
}

DualQuaternion DualQuaternion::operator+(const DualQuaternion &other) const {
	return DualQuaternion(real + other.real, dual + other.dual);
}

DualQuaternion DualQuaternion::operator-(const DualQuaternion &other) const {
	return DualQuaternion(real - other.real, dual - other.dual);
}

std::ostream& operator<<(std::ostream& ost, const DualQuaternion& toOut) {
	double realNorm, dualNorm;
	toOut.norm(&realNorm, &dualNorm);
	ost << toOut.real << ", " << toOut.dual << " (" << realNorm << ", " << dualNorm << ")";
	return ost;
}

void DualQuaternion::norm(double *realNorm, double *dualNorm) const {
	*realNorm = real.norm();
	*dualNorm = (real.w * dual.w + real.x * dual.x + real.y * dual.y + real.z * dual.z) / *realNorm;
}


DualQuaternion DualQuaternion::conjugate() const {
	return DualQuaternion(real.conjugate(), dual.conjugate());
}

DualQuaternion DualQuaternion::dualConjugate() const{
	return DualQuaternion(real, Quaternion(-dual.w, -dual.x, -dual.y, -dual.z));
}

void DualQuaternion::rotationPart(double *radians, double *x, double *y, double *z) {
	real.toAxisAngle(radians, x, y, z);
}

void DualQuaternion::rotationPart(double rotationMatrix[3][3]) {
	real.toRotationMatrix(rotationMatrix);
}

void DualQuaternion::translationPart(double *x, double *y, double *z) {
	Quaternion trans = dual * real.conjugate();
	*x = 2.0 * trans.x;
	*y = 2.0 * trans.y;
	*z = 2.0 * trans.z;
}

void DualQuaternion::translationPart(double vec[3]) {
	translationPart(vec, vec + 1, vec + 2);
}


}

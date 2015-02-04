#ifndef MGH_Quaternion_H
#define MGH_Quaternion_H

#include <iostream>

namespace LIB_NAMESPACE {
	class Quaternion {
		public:
			Quaternion();
			Quaternion(double initW, double initX, double initY, double initZ);

			Quaternion& operator=(const Quaternion &other);
			Quaternion operator+(const Quaternion &other) const;
			Quaternion operator-(const Quaternion &other) const;
			Quaternion operator*(const Quaternion &other) const;
			double norm() const;
			Quaternion normalize() const;
			Quaternion conjugate() const;

			void toAxisAngle(double *radians, double *x, double *y, double *z);
			static Quaternion fromAxisAngle(double radians, double x, double y, double z);
			void toRotationMatrix(double rotMatrix[3][3]);
			static Quaternion fromRotationMatrix(const double matrix[3][3]);


			friend std::ostream& operator<<(std::ostream& ost, const Quaternion& toOut);

			double w, x, y, z;

		private:
			static const double pi;
			static const double pix2;
	};


	class DualQuaternion {
		public:
			DualQuaternion();
			DualQuaternion(const Quaternion initReal, const Quaternion initDual);
			static DualQuaternion fromRotation(const double radians, const double x, const double y, const double z);
			static DualQuaternion fromRotationMatrix(const double matrix[3][3]);
			static DualQuaternion fromTranslation(const double x, const double y, const double z);

			void norm(double *realNorm, double *dualNorm) const;
			DualQuaternion conjugate() const;
			DualQuaternion dualConjugate() const;

			void rotationPart(double *radians, double *x, double *y, double *z);
			void rotationPart(double rotMatrix[3][3]);
			void translationPart(double *x, double *y, double *z);
			void translationPart(double transVec[3]);

			DualQuaternion& operator=(const DualQuaternion &other);
			DualQuaternion operator+(const DualQuaternion &other) const;
			DualQuaternion operator-(const DualQuaternion &other) const;
			DualQuaternion operator*(const DualQuaternion &other) const;
			friend std::ostream& operator<<(std::ostream& ost, const DualQuaternion& toOut);

			Quaternion real, dual;
	};
}

#endif

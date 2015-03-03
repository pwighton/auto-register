#ifndef MGH_AffineTransform_H
#define MGH_AffineTransform_H

#include <iostream>

namespace LIB_NAMESPACE {
	class AffineTransform {
		public:
			AffineTransform();
			AffineTransform(const double rotation[3][3], const double translation[3]);

			AffineTransform& operator=(const AffineTransform &other);
			AffineTransform operator*(const AffineTransform &other) const;
//			double* operator*(const double other[3]) const;

			AffineTransform inverse() const;

			void toAxisAngle(double *radians, double *x, double *y, double *z) const;
			static AffineTransform fromAxisAngleAndTranslation(double radians, double x, double y, double z, double tx, double ty, double tz);
			friend std::ostream& operator<<(std::ostream& ost, const AffineTransform& toOut);

			double matrix[4][4];

			void RotationPartZYXEulerAngles(double *eulerAngles) const;
			void RotationPartMatrix(double rot_matrix[3][3]) const;

			const static double identityRotation[3][3];
			const static double identityTranslation[3];
	};
}

#endif

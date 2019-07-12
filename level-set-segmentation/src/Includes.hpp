

#ifndef INCLUDES_HPP_
#define INCLUDES_HPP_

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <list>
#include <vector>
#include <stdlib.h>
#include <algorithm>
#include <inttypes.h>
#include <parallel/algorithm>
#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <list>
#include <ctime>
#include <time.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <inttypes.h>
#include <iomanip>
#include <locale>


#include <sys/stat.h>
#include <sys/time.h>
#include <chrono>
#include <opencv2/core/core.hpp>
#include <opencv2/core/mat.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/highgui/highgui_c.h>
#include <opencv2/ximgproc.hpp>
#include <opencv2/opencv.hpp>
#include <sys/stat.h>
#include <sys/time.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <chrono>

#include <Eigen/Dense>


#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include "opencv2/imgcodecs.hpp"
#include <opencv2/highgui.hpp>


using namespace cv;
using namespace std;
using namespace Eigen;

typedef uint64_t int_type_t;
typedef uint16_t uint_dist_type;
typedef uint16_t offset_int_type_t;

template<class T>
T FromString(const std::string& s)
{
	std::istringstream stream (s);
	T t;
	stream >> t;
	return t;
}

template<class T>
string ToString(T arg)
{
	std::ostringstream s;

	s << arg;

	return s.str();

}

enum DISTANCE {L1, L2, L2_approx, L_g, L2_induction };


inline void MultiplyMatrixVector(const vector< vector<double> >& M, int rows, int cols, const vector<double>& v, vector<double>& X){

	int c;
	// don't test to save time
	for (int r = 0; r < rows; r++){
		X[r] = 0;

		for (c = 0; c < cols; c++){
			X[r] += M[r][c]*v[c];

		}
	}


}

inline void MultiplySquareMatrixMatrix(const vector< vector<double> >& M0, const vector< vector<double> >& M1,
		int rows, vector< vector<double> >& R){

	int c;
	//double r;
	// don't test to save time
	for (int r = 0; r < rows; r++){


		for (c = 0; c < rows; c++){
			// row r in M0 * col c in M1.
			R[r][c] = 0;

			for (int inc = 0; inc < rows; inc++){
				R[r][c] += M0[r][inc]*M1[inc][c];
			}


		}
	}


}

inline void MultiplyMatricesWithSizes(const vector< vector<double> >& M0, const vector< vector<double> >& M1,
		int rowsA, int colsA, int colsB, vector< vector<double> >& R){

	int c;
	//double r;
	// don't test to save time
	for (int r = 0; r < rowsA; r++){


		for (c = 0; c < colsB; c++){
			// row r in M0 * col c in M1.
			R[r][c] = 0;

			for (int inc = 0; inc < colsA; inc++){
				R[r][c] += M0[r][inc]*M1[inc][c];
			}


		}
	}
}


inline void SubtractVectorFromVector(const vector<double>& A, const vector<double>& B, vector<double>& C, int rows){
	for (int r = 0; r < rows; r++){
		C[r] = A[r] - B[r];
	}
}

inline void AddVectorToVector(const vector<double>& A, const vector<double>& B, vector<double>& C, int rows){
	for (int r = 0; r < rows; r++){
		C[r] = A[r] + B[r];
	}
}

inline void AddMatrixToMatrix(const vector<vector< double> >& A, const vector<vector<double> >& B, vector<vector<double> >& C, int rows, int cols){
	for (int r = 0; r < rows; r++){
		for (int c = 0; c < cols; c++){
			C[r][c] = A[r][c] + B[r][c];
		}
	}
}

inline void MultiplyVectorByScalar(const vector<double>& A, const double scalar, vector<double>& B, int rows){
	for (int r = 0; r < rows; r++){
		B[r] = scalar*A[r];
	}

}

inline double SquaredDistance(const vector<double>& A, const vector<double>& B, int rows){
	double d = 0;
	for (int r = 0; r < rows; r++){
		d += (A[r] - B[r])*(A[r] - B[r]);
	}
	return d;
}

void PrintMatrix(vector< vector<double> >& p);



inline bool ProjectPointAndReturnIndex(const vector< vector<double> >& P,
		vector<double>& X, vector<double>& x, int rows, int cols, int_type_t& pixel_index){

	bool in = false;
	int r, c;
	MultiplyMatrixVector(P, 3, 4, X, x);


	x[0] /= x[2];  /// c
	x[1] /= x[2];  /// r
	x[2] = 1;


	c = round(x[0]);
	r = round(x[1]);

	if (c >= 0 && c < cols && r >= 0 && r < rows){
		in = true;
		pixel_index = round(x[1])*cols + round(x[0]);
	}	else {
		pixel_index = 0;
	}

	return in;
}

inline void NormalizePlane(vector<double>& p){

	double mag = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

	for (int i = 0; i < 4; i++){
		p[i] /= mag;
	}
}

inline void NormalizeVector(vector<double>& p){

	double mag = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

	if (fabs(mag) > 0.00000001){
	for (int i = 0; i < 3; i++){
		p[i] /= mag;
	}
	}
}

inline double MagnitudeVector(vector<double>& p){
	double mag = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
	return mag;
}

inline double DotProduct(const vector<double>& a, const vector<double>& b, int size){

	double r = 0;

	for (int i = 0; i < size; i++){
		r += a[i]*b[i];
	}
	return r;
}

inline void CrossProduct(const vector<double>& a, const vector<double>& b, vector<double>& c){

	c[0] = a[1]*b[2] - b[1]*a[2];
	c[1] = -a[0]*b[2] + b[0]*a[2];
	c[2] = a[0]*b[1] - b[0]*a[1];

}

inline void RayPlaneIntersection(const vector<double>& C, const vector<double>& V, const vector<double>& P, vector<double>& X ){

	double dp0 = DotProduct(C, P, 4);
	double dp1 = DotProduct(V, P, 4);

	if (dp1 != 0){
		double lambda = -dp0/dp1;

		X[0] = C[0] + lambda*V[0];
		X[1] = C[1] + lambda*V[1];
		X[2] = C[2] + lambda*V[2];

	}	else {

		X[0] = C[0];
		X[1] = C[1];
		X[2] = C[2];
		cout << "Slight error -- cannot find appropriate lambda b/c dot product is zero: " << endl;
		cout << "top " << dp0 << "   bottom " << dp1 << endl;
	}
}

inline bool RayPlaneIntersectionBool(const vector<double>& C, const vector<double>& V, const vector<double>& P, vector<double>& X ){

	double dp0 = DotProduct(C, P, 4);
	double dp1 = DotProduct(V, P, 4);

	if (dp1 != 0){
		double lambda = -dp0/dp1;

		X[0] = C[0] + lambda*V[0];
		X[1] = C[1] + lambda*V[1];
		X[2] = C[2] + lambda*V[2];

		if (lambda < 0){
			return false;
		}	else {
			return true;
		}
	}	else {

		X[0] = C[0];
		X[1] = C[1];
		X[2] = C[2];
		cout << "Slight error -- cannot find appropriate lambda b/c dot product is zero: " << endl;
		cout << "top " << dp0 << "   bottom " << dp1 << endl;

		return false;
	}
}




#endif /* INCLUDES_HPP_ */

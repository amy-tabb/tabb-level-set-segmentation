/*
 * LS0.hpp
 *
 *  Created on: Oct 17, 2017
 *      Author: atabb
 */

#ifndef IMAGEHISTOGRAM_HPP_
#define IMAGEHISTOGRAM_HPP_

//

//
//#include "SubSteps.hpp"
/// are all of these needed???
#include "Includes.hpp"
#include "ReconstructionStructure.hpp"
#include "DirectoryFunctions.hpp"
#include "DistanceTransforms.hpp"


class ImageHistogram{
public:
	vector<int_type_t> h;
	int_type_t n;
	double mu;
	double sigma_sq;

	ImageHistogram();

	~ImageHistogram();

	void AddToHistogram(int_type_t data_value);

	void Sum();

	void computeStatistics();

	void AddTo(const ImageHistogram& IH);

	void SubtractFrom(const ImageHistogram& IH);

	void Clear();

};

int IsNameInListLS(vector<string>& name_v, string test_name);

double ComputeEnergyWithHistograms(ImageHistogram& V1, ImageHistogram& V2, double v, int_type_t curve_length, int_type_t number_aged_out );


int_type_t indexRC(int_type_t r, int_type_t c, int rows);

void XYZFromIndex(int_type_t index, int_type_t& x, int_type_t& y, int_type_t& z, int_type_t xsize, int_type_t ysize);

int_type_t ReturnIndexFromXYZ(int_type_t x, int_type_t y, int_type_t z, int_type_t xsize, int_type_t ysize);

double ComputeEnergyWithHistograms(vector<ImageHistogram>& V1, vector<ImageHistogram>& V2, double v, int_type_t curve_length, int aged_out);




#endif /* IMAGEHISTOGRAM_HPP_ */

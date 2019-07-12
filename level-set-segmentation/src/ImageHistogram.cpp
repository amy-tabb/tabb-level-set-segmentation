/*
 * LS0.cpp
 *
 *  Created on: Oct 17, 2017
 *      Author: atabb
 */




#include "ImageHistogram.hpp"

ImageHistogram::ImageHistogram(){

	h.resize(256, 0);
}

ImageHistogram::~ImageHistogram(){

}

void ImageHistogram::AddToHistogram(int_type_t data_value){
	if (data_value < 256){
		h[data_value]++;
	}	else {
		cout << "Data value out of bounds " << data_value << endl;
		exit(1);
	}
}

void ImageHistogram::Sum(){
	n = 0;
	for (int i = 0; i < 256; i++){
		n = n + h[i];
	}
}

void ImageHistogram::computeStatistics(){
	Sum();  /// make sure this is up to date.
	cout << "Compute statistics :, n " << n << endl;
	mu = 0;

	for (int i = 0; i < 256; i++){
		mu = mu + double(h[i])*double(i)/255.0;
	}
	if (n > 0){
		mu = mu/double(n);


		sigma_sq = 0;
		for (int i = 0; i < 256; i++){
			sigma_sq = sigma_sq + double(h[i])*(pow(double(i)/255.0 - mu, 2))/double(n);
		}
	}	else {
		sigma_sq = 0;
	}
}

int IsNameInListLS(vector<string>& name_v, string test_name){

	for (int j = 0, jn = name_v.size(); j < jn; j++){
		if (strcmp(test_name.c_str(), name_v[j].c_str()) == 0){
			return j;
		}
	}

	return -1;

}

int_type_t indexRC(int_type_t r, int_type_t c, int rows){
	return c*rows + r;
}



double ComputeEnergyWithHistograms(ImageHistogram& V1, ImageHistogram& V2, double v, int_type_t curve_length, int_type_t number_aged_out ){



	//int_type_t rs_index = 0;
	double current_energy= 0;
	double pdf_1 = 0;
	double pdf_2 = 0;

//	double denominator_class1 = 1.0/(sqrt(2.0*3.14*V1.sigma_sq));
//	double denominator_class2 = 1.0/(sqrt(2.0*3.14*V2.sigma_sq));

	double class1_e = 0;
	double class2_e = 0;

	cout << "Stats " << V1.mu << " " << V1.sigma_sq << ", " << V2. mu << ", " << V2.sigma_sq << endl;

	double id;
	for (int i = 0; i < 256; i++){

		id = double(i)/255.0;
		pdf_1 = exp(-pow(id - V1.mu, 2)/(2.0*V1.sigma_sq))/(sqrt(2.0*3.14*V1.sigma_sq));
		// how many are there?
		if (pdf_1 > 0){
			class1_e += double(V1.h[i])*log10(pdf_1);
		}


		double diff2 = id - V2.mu;
		double exponent = -(diff2*diff2)/(2*V2.sigma_sq);
		double numer = exp(exponent);
		if (numer > 0){
			pdf_2 = numer/(sqrt(2.0*3.14*V2.sigma_sq));
			class2_e += double(V2.h[i])*log10(pdf_2);
		}	else {
			pdf_2 = 0;
		}
	}


	current_energy = -class1_e + -class2_e + v*double(curve_length + number_aged_out);
	return current_energy;
}



void ImageHistogram::AddTo(const ImageHistogram& IH){
	for (int i = 0; i < 256; i++){
		h[i] += IH.h[i];
	}

}

void ImageHistogram::SubtractFrom(const ImageHistogram& IH){
	for (int i = 0; i < 256; i++){
		if (IH.h[i] > h[i]){
			cout << "ERROR about to have an overflow in a histogram, subtract " << endl;
			cout << "Main " << h[i] << " to subtract  " << IH.h[i] << endl;
			cout << " on value " << i << endl;
			exit(1);
		}
		h[i] -= IH.h[i];
	}

}

void ImageHistogram::Clear(){
	for (int i = 0; i < 256; i++){
		h[i] = 0;
	}

	mu = 0;
	sigma_sq = 0;
	n = 0;
}


double ComputeEnergyWithHistograms(vector<ImageHistogram>& V1, vector<ImageHistogram>& V2, double v, int_type_t curve_length, int aged_out){

	double summed_energy = 0;


	for (int h = 0, in = V1.size(); h < in; h++){

		double current_energy= 0;

		double denominator_class1 = 1.0;
		double denominator_class2 = 1.0;

		double pdf_1 = 0;
		double pdf_2 = 0;

		double class1_e = 0;
		double class2_e = 0;



		//cout << "Stats " << V1.mu << " " << V1.sigma_sq << ", " << V2. mu << ", " << V2.sigma_sq << endl;


		double id;

		if (V1[h].n > 0){
			denominator_class1 = (sqrt(2.0*3.14*V1[h].sigma_sq));

			for (int i = 0; i < 256; i++){

				id = double(i)/255.0;
				//pdf_1 = exp(-pow(id - V1[h].mu, 2)/(2.0*V1[h].sigma_sq))/(sqrt(2.0*3.14*V1[h].sigma_sq));
				pdf_1 = exp(-pow(id - V1[h].mu, 2)/(2.0*V1[h].sigma_sq))/denominator_class1;
				// how many are there?
				if (pdf_1 > 0){
					class1_e += double(V1[h].h[i])*log10(pdf_1);
				}

			}
		}


		if (V2[h].n > 0){
			denominator_class2 = (sqrt(2.0*3.14*V2[h].sigma_sq));

			for (int i = 0; i < 256; i++){
				//pdf_2 = exp(-pow(id - V2[h].mu, 2)/(2.0*V2[h].sigma_sq))/(sqrt(2.0*3.14*V2[h].sigma_sq));
				pdf_2 = exp(-pow(id - V2[h].mu, 2)/(2.0*V2[h].sigma_sq))/denominator_class2;
				// how many are there?
				if (pdf_2 > 0){
					class2_e += double(V2[h].h[i])*log10(pdf_2);
				}
			}
		}

		current_energy = -class1_e + -class2_e + v*double(curve_length + aged_out);

		summed_energy += current_energy;

	}

	return summed_energy;
}
// z dominant version here.
void XYZFromIndex(int_type_t index, int_type_t& x, int_type_t& y, int_type_t& z, int_type_t xsize, int_type_t ysize){
	z = index/((ysize*xsize));
	int_type_t temp = index % ((ysize*xsize));
	x = temp / ysize;
	y = temp % ysize;

}

// z dominant version here.
int_type_t ReturnIndexFromXYZ(int_type_t x, int_type_t y, int_type_t z, int_type_t xsize, int_type_t ysize){
	return z*(ysize*xsize) + x*ysize + y;
}




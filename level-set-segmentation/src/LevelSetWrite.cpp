/*
 * LevelSetWrite.cpp
 *
 *  Created on: May 3, 2019
 *      Author: atabb
 */

#include "LevelSets.hpp"


void WriteResultOnImages(string iter_write_dir, int_type_t number_images, vector<Mat>& images, vector<bool*>& class_flag){


	string command = "mkdir " + iter_write_dir;
	string filename;
	Mat cloned_image;
	Vec3b intensity_bgr;

	system(command.c_str());

	int rows = images[0].rows;
	int cols = images[0].cols;


#pragma omp parallel for private(filename, cloned_image, intensity_bgr)
	for (int_type_t i = 0; i < number_images; i++){
		//thread_id = omp_get_thread_num();
		cvtColor(images[i], cloned_image, cv::COLOR_GRAY2BGR);

		for (int r = 0; r < rows; r++){
			for (int c = 0; c < cols; c++){
				if (class_flag[i][c*rows + r] == false){
					intensity_bgr = cloned_image.at<Vec3b>(r, c);
					intensity_bgr[0] = 255;
					cloned_image.at<Vec3b>(r, c) = intensity_bgr;
				}
			}

		}

		filename = ToString<int>(i);

		while (filename.size() < 5){
			filename = "0" + filename;
		}

		filename = iter_write_dir + "/" + filename + ".png";
		imwrite(filename.c_str(), cloned_image);


	}



}

void WriteThresholdsOnImages(string iter_write_dir, int_type_t number_images, vector<Mat>& images, vector<bool*>& class_flag, int lower_threshold, int upper_threshold,
		int background_threshold){
	string command = "mkdir " + iter_write_dir;
	int thread_id;
	Vec3b intensity_bgr;
	int image_value;
	string filename;

	system(command.c_str());
	Mat cloned_image;

	int rows = images[0].rows;
	int cols = images[0].cols;


#pragma omp parallel for private(filename, cloned_image, image_value, intensity_bgr)
	for (int_type_t i = 0; i < number_images; i++){
		//thread_id = omp_get_thread_num();
		cvtColor(images[i], cloned_image, cv::COLOR_GRAY2BGR);

		for (int r = 0; r < rows; r++){
			for (int c = 0; c < cols; c++){

				//image_value = cloned_image.data[3*(r*cols + c) + 0] ;
				intensity_bgr = cloned_image.at<Vec3b>(r, c);
				image_value = intensity_bgr[0];

				if (!(image_value > lower_threshold && image_value < upper_threshold)){

					//cloned_image.data[3*(r*cols + c) + 0] = 255;
					intensity_bgr[0] = 255;


				}	else {
					//cloned_image.data[3*(r*cols + c) + 1] = 255;
					intensity_bgr[1] = 255;
				}
				cloned_image.at<Vec3b>(r, c) = intensity_bgr;

			}

		}

		filename = ToString<int>(i);

		while (filename.size() < 5){
			filename = "0" + filename;
		}

		filename = iter_write_dir + "/" + filename + ".png";
		imwrite(filename.c_str(), cloned_image);


	}

}


void WriteBWResultImages(string iter_write_dir, int_type_t number_images, vector<Mat>& bw_ims, vector<bool*>& class_flag){

	string command = "mkdir " + iter_write_dir;
	int thread_id;

	system(command.c_str());

	int rows = bw_ims[0].rows;
	int cols = bw_ims[0].cols;
	string filename;

#pragma omp parallel for private(thread_id, filename)
	for (int_type_t i = 0; i < number_images; i++){
		thread_id = omp_get_thread_num();
		bw_ims[thread_id].setTo(0);

		for (int r = 0; r < rows; r++){
			for (int c = 0; c < cols; c++){
				if (class_flag[i][c*rows + r] == false){
					bw_ims[thread_id].at<uchar>(r, c) = 255;
				}
			}

		}

		filename = ToString<int>(i);

		while (filename.size() < 5){
			filename = "0" + filename;
		}

		filename = iter_write_dir + "/" + filename + ".png";
		imwrite(filename.c_str(), bw_ims[thread_id]);

	}

}



void WriteBWResultImagesConnectedComponents(string iter_write_dir, int_type_t number_images, vector<Mat>& bw_ims, vector<bool*>& class_flag, bool* cc_we_want,
		vector<uint32_t*>& cc_map, int index_biggest_cc ){

	string command = "mkdir " + iter_write_dir;
	int thread_id;

	system(command.c_str());

	int rows = bw_ims[0].rows;
	int cols = bw_ims[0].cols;
	string filename;

	uint32_t cc_value;

	if (index_biggest_cc == -1){



		#pragma omp parallel for private(thread_id, filename, cc_value)
		for (int i = 0; i < number_images; i++){
			thread_id = omp_get_thread_num();
			bw_ims[thread_id].setTo(0);
			for (int r = 0; r < rows; r++){
				for (int c = 0; c < cols; c++){

					cc_value  = cc_map[i][c*rows + r];

					if (cc_we_want[cc_value] == true){
						bw_ims[thread_id].at<uchar>(r, c) = 255;
					}
				}

			}

			filename = ToString<int>(i);

			while (filename.size() < 5){
				filename = "0" + filename;
			}

			filename = iter_write_dir + "/" + filename + ".png";
			imwrite(filename.c_str(), bw_ims[thread_id]);

		}
	}	else {

		#pragma omp parallel for private(thread_id, filename, cc_value)
				for (int i = 0; i < number_images; i++){
					thread_id = omp_get_thread_num();
					bw_ims[thread_id].setTo(0);
					if (index_biggest_cc >= 0){
						for (int r = 0; r < rows; r++){
							for (int c = 0; c < cols; c++){

								cc_value  = cc_map[i][c*rows + r];

								if (cc_value == int_type_t(index_biggest_cc)){
									//bw_ims[thread_id].data[r*cols + c] = 255;
									bw_ims[thread_id].at<uchar>(r, c) = 255;
								}
							}

						}
					}

					filename = ToString<int>(i);

					while (filename.size() < 5){
						filename = "0" + filename;
					}

					filename = iter_write_dir + "/" + filename + ".png";
					imwrite(filename.c_str(), bw_ims[thread_id]);

				}
	}

}


void WriteGridImages(string iter_write_dir, int_type_t number_images, vector<Mat>& bw_ims, vector<bool*>& grid_active, int grid_resolution,
		int_type_t grid_xsize, int_type_t grid_ysize, int_type_t grid_zsize){


	int_type_t gx, gy, gz;
	string command = "mkdir " + iter_write_dir;
	int thread_id;
	string filename;

	system(command.c_str());

	int rows = bw_ims[0].rows;
	int cols = bw_ims[0].cols;


#pragma omp parallel for private(thread_id, filename, gx, gy, gz)
	for (gz=  0; gz < grid_zsize; gz++){
		thread_id = omp_get_thread_num();
		bw_ims[thread_id].setTo(0);
		for (gx=  0; gx < grid_xsize; gx++){
			for (gy=  0; gy < grid_ysize; gy++){
				if (grid_active[gz][gx*grid_ysize + gy] == true){
					for (int_type_t r = gy*grid_resolution; r < int_type_t(rows) && r < (gy+1)*grid_resolution; r++){
						for (int_type_t c = gx*grid_resolution; c < int_type_t(cols) && c < (gx+1)*grid_resolution; c++){

							//bw_ims[thread_id].data[r*cols + c] = 255;
							bw_ims[thread_id].at<uchar>(r, c) = 255;
						}
					}
				}

			}
		}

		for (int_type_t i = gz*grid_resolution; i < int_type_t(number_images) && i < (gz+1)*grid_resolution; i++){
			filename = ToString<int>(i);

			while (filename.size() < 5){
				filename = "0" + filename;
			}

			filename = iter_write_dir + "/" + filename + ".png";
			imwrite(filename.c_str(), bw_ims[thread_id]);

		}

	}


}

void WriteDistanceImages(string iter_write_dir, int_type_t number_images, vector<Mat>& bw_ims, vector<uint_dist_type*>& dt_xyz,
		vector<bool*>& class_flag, vector<int8_t*>& contour_age, int contour_age_out, uint_dist_type band_size_sq){

	string command = "mkdir " + iter_write_dir;
	int thread_id;
	string filename;

	system(command.c_str());


	int rows = bw_ims[0].rows;
	int cols = bw_ims[0].cols;

	// could be parallelized.
#pragma omp parallel for private(thread_id, filename)
	for (int i = 0; i < number_images; i++){
		thread_id = omp_get_thread_num();
		bw_ims[thread_id].setTo(0);
		for (int r = 0; r < rows; r++){
			for (int c = 0; c < cols; c++){

				if (class_flag[i][c*rows + r] == false){
					if (dt_xyz[i][c*rows + r] == 0){ /// only on the active contour. -- want to see the active contour.

						//bw_ims[thread_id].data[r*cols + c] = 250;
						bw_ims[thread_id].at<uchar>(r, c) = 250;

						if (contour_age[i][c*rows + r] >= contour_age_out){ /// since age out is 1, this shjould be evceryone.
							//bw_ims[thread_id].data[r*cols + c] = 200;
							bw_ims[thread_id].at<uchar>(r, c) = 200;
						}
					}	else {
						if (contour_age[i][c*rows + r] >= contour_age_out){
							//bw_ims[thread_id].data[r*cols + c] = 150;
							bw_ims[thread_id].at<uchar>(r, c) = 150;
						}	else {
							//bw_ims[thread_id].data[r*cols + c] = 100;
							bw_ims[thread_id].at<uchar>(r, c) = 100;
						}
					}
				}	else {

					if (dt_xyz[i][c*rows + r] <= band_size_sq){
						//bw_ims[thread_id].data[r*cols + c] = 50;
						bw_ims[thread_id].at<uchar>(r, c) = 50;
					}
				}

			}

		}

		filename = ToString<int>(i);

		while (filename.size() < 5){
			filename = "0" + filename;
		}

		filename = iter_write_dir + "/" + filename + ".png";
		imwrite(filename.c_str(), bw_ims[thread_id]);

	}
}

/*
 * LevelSets.hpp
 *
 *  Created on: May 2, 2019
 *      Author: atabb
 */

#ifndef LEVELSETS_HPP_
#define LEVELSETS_HPP_

#include "ImageHistogram.hpp"

void LevelSetMainFunction(string test_image_dir, string result_image_dir, int lower_threshold, int upper_threshold,
		int background_lower_threshold, double v, int grid_resolution, uint_dist_type band_size,
		bool mark_all, int contour_age_out, int min_active_contour_size, int max_number_threads, bool write_on_image, bool write_initial );

void CurvatureAndDTParallelOneAllocLayers(vector<bool*>& class_flag, vector<uint_dist_type*>& dt_xyz, uint_dist_type band_size,
		bool two_d_images,  int_type_t xsize, int_type_t ysize, int_type_t zsize, vector<Mat>& images,
		ImageHistogram& V1, ImageHistogram& V2, ImageHistogram& V12, ImageHistogram& V21, vector< float* >& current_layers,
		float v, float speed, int lower_threshold, int upper_threshold,
		vector<bool*>& active_grid, int_type_t grid_xsize, int_type_t grid_ysize, int_type_t grid_resolution );

float computeCurvatureFloat(vector<float*>& phi, int r, int c, int i, bool two_d_curvature, int rows, float& maxgrad, float& mingrad);

void WriteResultOnImages(string iter_write_dir, int_type_t number_images, vector<Mat>& images, vector<bool*>& class_flag);

void WriteThresholdsOnImages(string iter_write_dir, int_type_t number_images, vector<Mat>& images, vector<bool*>& class_flag, int lower_threshold, int upper_threshold,
		int background_lower_threshold);

void WriteBWResultImages(string iter_write_dir, int_type_t number_images, vector<Mat>& bw_ims, vector<bool*>& class_flag);

void WriteGridImages(string iter_write_dir, int_type_t number_images, vector<Mat>& bw_ims, vector<bool*>& grid_active, int grid_resolution,
		int_type_t grid_xsize, int_type_t grid_ysize, int_type_t grid_zsize);


void WriteDistanceImages(string iter_write_dir, int_type_t number_images, vector<Mat>& bw_ims, vector<uint_dist_type*>& dt_xyz,
		vector<bool*>& class_flag, vector<int8_t*>& contour_age, int contour_age_out, uint_dist_type band_size_sq);

void WriteBWResultImagesConnectedComponents(string iter_write_dir, int_type_t number_images, vector<Mat>& bw_ims, vector<bool*>& class_flag, bool* cc_we_want,
		vector<uint32_t*>& cc_map, int index_biggest_cc );

pair<int_type_t, int_type_t> UpdateEdgeListDiffStructure(vector<bool*>& class_flag, vector<uint_dist_type*>& dt_xyz,
		vector<int8_t*>& contour_age, vector<vector< int_type_t*> >& d_gamma_indices_per_thread,
		vector<int_type_t>& count_d_gamma_per_thread, int number_threads, int_type_t xsize, int_type_t ysize,
		int_type_t zsize,
		bool use_band_size, uint_dist_type band_size_sq, int8_t contour_age_threshold);

void UpdateGridWithNewEdges(vector<vector< int_type_t*> >& d_gamma_indices_per_thread,
		vector<int_type_t>& count_d_gamma_per_thread, vector<bool*>& grid, int_type_t xsize, int_type_t ysize,
		int_type_t grid_xsize, int_type_t grid_ysize, int_type_t grid_zsize, int_type_t grid_resolution, int_type_t number_threads);

void ReturnLimits(int_type_t gx, int_type_t gy, int_type_t gz, int_type_t mx, int_type_t my, int_type_t mz,
		int_type_t& minx, int_type_t& miny, int_type_t& minz,
		int_type_t& maxx, int_type_t& maxy, int_type_t& maxz);

void ConnectedComponents(vector<uint32_t*>& cc_map, vector<int_type_t>& counts_per, int xsize, int ysize, int zsize);

int_type_t InitiateSearch(vector<uint32_t*>& cc_map, int_type_t xi, int_type_t yi, int_type_t zi, int_type_t xsize, int_type_t ysize, int_type_t zsize, uint32_t cc_number);






#endif /* LEVELSETS_HPP_ */

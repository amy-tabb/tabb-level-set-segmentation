
#ifndef DISTANCETRANSFORMS_HPP_
#define DISTANCETRANSFORMS_HPP_


#include "Includes.hpp"
#include "ReconstructionStructure.hpp"


void DistanceTransformSqdParallelTubeXIV(vector<vector< int_type_t*> >& d_gamma_indices_per_thread,
		vector<int_type_t>& count_d_gamma_per_thread, vector<uint_dist_type*>& dt_xyz, uint_dist_type big_number,
		uint_dist_type band_increment, int_type_t xsize, int_type_t ysize, int_type_t zsize, vector<bool*>& grid, int_type_t grid_xsize,
		int_type_t grid_ysize, int_type_t grid_zsize, int_type_t grid_resolution, int_type_t number_threads);

void ReturnXYZIndicesFromIndex(int_type_t index, int_type_t& x, int_type_t& y, int_type_t& z, int_type_t xsize, int_type_t ysize);


#endif /* DISTANCETRANSFORMS_HPP_ */

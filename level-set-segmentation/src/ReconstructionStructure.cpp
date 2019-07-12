#include "ReconstructionStructure.hpp"


ReconstructionStructure::ReconstructionStructure(){
	number_voxels_grid = 0;
}

ReconstructionStructure::~ReconstructionStructure(){
}

void ReconstructionStructure::CreateGridMinimal(vector< vector<double> >& boundingbox, double division){

	inter_voxel_distance = division;
	BB = boundingbox;
	initial_offset.resize(3, 0);
	number_voxels_per_dim.resize(3, 0);

	initial_offset[0] = BB[0][0] + inter_voxel_distance/2.0;
	initial_offset[1] = BB[0][1] + inter_voxel_distance/2.0;
	initial_offset[2] = BB[0][2] + inter_voxel_distance/2.0;

	for (int i = 0; i < 3; i++){
		number_voxels_per_dim[i] = floor((BB[1][i] - BB[0][i])/inter_voxel_distance);
		cout << "Number voxels per " << number_voxels_per_dim[i] << endl;
	}

	// configuration grid is only for the voxels.
	number_voxels_grid = number_voxels_per_dim[0]*number_voxels_per_dim[1]*number_voxels_per_dim[2];

	// later - make faster, make an array for ALL of these things ...
	configuration_grid.resize(number_voxels_grid, true);
}



int_type_t ReconstructionStructure::ReturnIndexFromXYZIndices(int_type_t x, int_type_t y, int_type_t z){
	if (x*(number_voxels_per_dim[1])*(number_voxels_per_dim[2]) + y*(number_voxels_per_dim[2]) + z >= number_voxels_grid){
		cout << "ERROR ON size " <<x*(number_voxels_per_dim[1])*(number_voxels_per_dim[2]) + y*(number_voxels_per_dim[2]) + z << ", n voxesl " << number_voxels_grid << endl;
		cout << "x , y, z " << x << ", " << y << ",  " << z << endl;
		exit(1);
	}

	return x*(number_voxels_per_dim[1])*(number_voxels_per_dim[2]) + y*(number_voxels_per_dim[2]) + z;

}

bool ReconstructionStructure::SizeCheck(int_type_t x, int_type_t y, int_type_t z){
	if (x >= number_voxels_per_dim[0]){
		cout << "x is too big " << x << " out of possible " << number_voxels_per_dim[0] << endl;
		return false;
	}

	if (y >= number_voxels_per_dim[1]){
			cout << "y is too big " << y << " out of possible " << number_voxels_per_dim[1] << endl;
			return false;
		}

	if (z >= number_voxels_per_dim[2]){
			cout << "z is too big " << z << " out of possible " << number_voxels_per_dim[2] << endl;
			return false;
		}

	return true;



}

void ReconstructionStructure::ReturnXYZIndicesFromIndex(int_type_t voxel_index, int_type_t& x, int_type_t& y, int_type_t& z){

	int_type_t temp_index = voxel_index;
	x = temp_index/(number_voxels_per_dim[1]*number_voxels_per_dim[2]);

	temp_index -= x*(number_voxels_per_dim[1])*(number_voxels_per_dim[2]);
	y = temp_index/(number_voxels_per_dim[2]);

	temp_index -= y*(number_voxels_per_dim[2]);

	z = temp_index;

	if (x*(number_voxels_per_dim[1])*(number_voxels_per_dim[2]) + y*(number_voxels_per_dim[2]) + z != voxel_index){
		cout << "ERROR on vox computation! " << endl << endl;
		exit(1);
	}
}

int ReconstructionStructure::Return26ConnectedNeighbors(int_type_t start_voxel, int_type_t* member_array, int_type_t* distance){

	int number_neighbors = 0;
	int_type_t nindex;

	int_type_t x_index, y_index, z_index;

	ReturnXYZIndicesFromIndex(start_voxel, x_index, y_index, z_index);


	int_type_t minx, maxx, miny, maxy, minz, maxz;

	if (x_index == 0){
		minx = 0;
	}	else {
		minx = x_index - 1;
	}

	if (y_index == 0){
		miny = 0;
	}	else {
		miny = y_index - 1;
	}

	if (z_index == 0){
		minz = 0;
	}	else {
		minz = z_index - 1;
	}

	if (z_index == number_voxels_per_dim[2] - 1){
		maxz = z_index;
	}	else {
		maxz = z_index + 1;
	}

	if (y_index == number_voxels_per_dim[1] - 1){
		maxy = y_index;
	}	else {
		maxy = y_index + 1;
	}

	if (x_index == number_voxels_per_dim[0] - 1){
		maxx = x_index;
	}	else {
		maxx = x_index + 1;
	}

	for (int_type_t x0 = minx; x0 <= maxx; x0++){
		for (int_type_t y0 = miny; y0 <= maxy; y0++){
			for (int_type_t z0 = minz; z0 <= maxz; z0++){
				nindex = ReturnIndexFromXYZIndices(x0, y0, z0);
				if (nindex != start_voxel){
					member_array[number_neighbors] = nindex;
					distance[number_neighbors] = (x0 != x_index) + (y0 != y_index) + (z0 != z_index);
					number_neighbors++;
				}
			}
		}
	}


	return number_neighbors;


}



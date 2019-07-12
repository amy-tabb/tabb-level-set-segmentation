#ifndef RECONSTRUCTIONSTRUCTURE_HPP_
#define RECONSTRUCTIONSTRUCTURE_HPP_
#include "Includes.hpp"


class ReconstructionStructure{
public:
	vector<vector<double > > all_planes;

	vector<bool> configuration_grid;

	vector<double> initial_offset;
	double inter_voxel_distance;
	vector<int_type_t> number_voxels_per_dim;
	int_type_t number_voxels_grid;

	vector< vector<double> > BB;


	ReconstructionStructure();

	~ReconstructionStructure();

	int_type_t ReturnIndexFromXYZIndices(int_type_t x, int_type_t y, int_type_t z);

	void WritePlyFile(string outfile,
			vector<int_type_t>& subdims,
			vector<double>& points,
			vector<int_type_t>& faces,
			vector<int_type_t>& color);

	void ReturnXYZIndicesFromIndex(int_type_t voxel_index, int_type_t& x, int_type_t& y, int_type_t& z);

	void CreateGridMinimal(vector< vector<double> >& boundingbox, double division);


	void GenerateAndWriteSurfaceInPlyFormat(string outdir, int iteration,
			string prefix, vector<int>& cs, bool default_string);

	void GenerateAndWriteSurfaceInPlyFormat(string outdir, int iteration, string prefix = "iter", int* cs=NULL, bool default_string = true);

	bool SizeCheck(int_type_t x, int_type_t y, int_type_t z);

	int Return26ConnectedNeighbors(int_type_t start_voxel, int_type_t* member_array, int_type_t* distance);

};


template<class T>
std::string FormatWithCommas(T value)
{
	int_type_t uvalue = value;
	bool negative = false;
	if (value < 0){
		negative = true;
		uvalue = -value;
	}


	string s;
	  int cnt = 0;
	  do
	  {
	    s.insert(0, 1, char('0' + uvalue % 10));
	    uvalue /= 10;
	    if (++cnt == 3 && uvalue)
	    {
	      s.insert(0, 1, ',');
	      cnt = 0;
	    }
	  } while (uvalue);

	  if (negative){
		  s = "-" + s;
	  }
	  return s;
}

void PrintVector(vector<double>& p);

bool IsLegalPoint(ReconstructionStructure& RS, int source_cam, int within_cam_index, pair<int, int>cam_label_pair);



#endif /* RECONSTRUCTIONSTRUCTURE_HPP_ */

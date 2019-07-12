/*
 * LevelSetsConnectedComponents.cpp
 *
 *  Created on: May 3, 2019
 *      Author: atabb
 */


#include "LevelSets.hpp"

void ConnectedComponents(vector<uint32_t*>& cc_map, vector<int_type_t>& counts_per, int xsize, int ysize, int zsize){

	int_type_t current_cc = 2;
	counts_per.resize(2, 0);
	int_type_t count = 0;

	for (int zi = 0; zi < zsize; zi++){
		/// compute statistics by image, count how many per class per image.
		for (int xi = 0; xi < xsize; xi++){
			for (int yi = 0; yi < ysize; yi++){

				if (cc_map[zi][xi*ysize + yi] == 0){
					// call other function.
					count = InitiateSearch(cc_map, xi, yi, zi, xsize, ysize, zsize, current_cc);

					counts_per.push_back(count);
					current_cc++;

				}

			}
		}
	}
}

int_type_t InitiateSearch(vector<uint32_t*>& cc_map, int_type_t xi, int_type_t yi, int_type_t zi, int_type_t xsize, int_type_t ysize, int_type_t zsize, uint32_t cc_number){

	int_type_t count = 0;
	vector<int_type_t> items_in_cc;
	int_type_t current_index, neighbor_index;

	current_index = ReturnIndexFromXYZ(xi, yi, zi, xsize, ysize);

	items_in_cc.push_back(current_index);

	cc_map[zi][xi*ysize + yi] = cc_number;

	int_type_t zs, ze, xs, xe, ys, ye;


	for (  ; count < items_in_cc.size(); count++){
		current_index = items_in_cc[count];
		XYZFromIndex(current_index, xi, yi, zi, xsize, ysize);

		zi == 0 ? zs = 0 : zs = zi - 1;
		yi == 0 ? ys = 0 : ys = yi - 1;
		xi == 0 ? xs = 0 : xs = xi - 1;

		zi == zsize - 1 ? ze = zsize - 1 : ze = zi + 1;
		yi == ysize - 1 ? ye = ysize - 1 : ye = yi + 1;
		xi == xsize - 1 ? xe = xsize - 1 : xe = xi + 1;



		for (int_type_t z0 = zs; z0 <= ze; z0++){
			for (int_type_t x0 = xs; x0 <= xe; x0++){
				for (int_type_t y0 = ys; y0 <= ye; y0++){
					if (cc_map[z0][x0*ysize + y0] == false){  // occupied, but not yet labeled.
						cc_map[z0][x0*ysize + y0] = cc_number;
						neighbor_index = ReturnIndexFromXYZ(x0, y0, z0, xsize, ysize);
						items_in_cc.push_back(neighbor_index);
					}
				}
			}
		}

	}

	return count;

}

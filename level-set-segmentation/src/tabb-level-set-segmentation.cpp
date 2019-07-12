//This README file is to accompany code for robot-world, hand-eye calibration, produced by Amy Tabb as companion to a paper:
//	Solving the Robot-World Hand-Eye(s) Calibration Problem with Iterative Methods
//
//```
//@inproceedings{tabb_segmenting_2018,
//	title = {Segmenting {Root} {Systems} in {X}-{Ray} {Computed} {Tomography} {Images} {Using} {Level} {Sets}},
//	doi = {10.1109/WACV.2018.00070},
//	abstract = {The segmentation of plant roots from soil and other growing media in X-ray computed tomography images is needed to effectively study the root system architecture without excavation. However, segmentation is a challenging problem in this context because the root and non-root regions share similar features. In this paper, we describe a method based on level sets and specifically adapted for this segmentation problem. In particular, we deal with the issues of using a level sets approach on large image volumes for root segmentation, and track active regions of the front using an occupancy grid. This method allows for straightforward modifications to a narrow-band algorithm such that excessive forward and backward movements of the front can be avoided, distance map computations in a narrow band context can be done in linear time through modification of Meijster et al.'s distance transform algorithm, and regions of the image volume are iteratively used to estimate distributions for root versus non-root classes. Results are shown of three plant species of different maturity levels, grown in three different media. Our method compares favorably to a state-of-the-art method for root segmentation in X-ray CT image volumes.},
//	booktitle = {2018 {IEEE} {Winter} {Conference} on {Applications} of {Computer} {Vision} ({WACV})},
//	author = {Tabb, A. and Duncan, K. E. and Topp, C. N.},
//	month = mar,
//	year = {2018},
//	keywords = {Shape, Image segmentation, Computed tomography, Level set, Media, Soil, X-ray imaging},
//	pages = {586--595},
//}
//```
// license: MIT, no warranties expressed or implied.
//MIT License
//
//Copyright (c) 2019 Amy Tabb
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

// For updates, please check the Github repository: https://github.com/amy-tabb/tabb-level-set-segmentation
// July 2019

#include <sys/stat.h>
#include "Includes.hpp"
#include "ReconstructionStructure.hpp"
#include "DirectoryFunctions.hpp"
#include "DistanceTransforms.hpp"
#include "LevelSets.hpp"
#include "ConvertRepresentation.hpp"
#include <getopt.h>

void EnsureDirHasTrailingBackslash(string& write_directory){
	int n_letters = write_directory.size();
	bool eval =  (write_directory[n_letters - 1] == '/');
	cout << "Last character compare " << write_directory << " " <<  eval << endl;
	if (eval == false){
		write_directory = write_directory + "/";
	}

}



int main (int argc, char **argv)
{
	int do_preprocessing = false;
	int do_segmentation = false;
	int print_help = false;

	string input_directory;
	string output_directory;
	int lower_threshold = 0;
	int upper_threshold = 255;
	int background_threshold = 0;
	double nu = 1.0;
	int grid_resolution = 10;
	int band_size = 10;
	bool disregard_histogram_grid = false;
	int t = 1;
	int min_active_contour_size = 100;
	int max_number_threads = omp_get_max_threads();
	bool write_on_image = false;
	bool write_initial = false;


	int opt_argument;

	while (1)
	{
		static struct option long_options[] =
		{
				/* These options set a flag. */
				{"preprocessing", no_argument,       &do_preprocessing, 1},
				{"segmentation",   no_argument,       &do_segmentation, 1},
				{"help",   no_argument,       &print_help, 1},
				/* These options don’t set a flag.
             We distinguish them by their indices. */
				{"input",   required_argument, 0, 'a'},
				{"output",  required_argument, 0, 'b'},
				{"lowerthresh",  required_argument, 0, 'c'},
				{"upperthresh",  required_argument, 0, 'd'},
				{"background",    required_argument, 0, 'e'},
				{"nu",  required_argument, 0, 'f'},
				{"grid-resolution",    required_argument, 0, 'g'},
				{"band-size",    required_argument, 0, 'h'},
				{"disregard-histogram-grid",    required_argument, 0, 'i'},
				{"t",    required_argument, 0, 'j'},
				{"min-contour-size", required_argument, 0, 'k'},
				{"max-threads", required_argument, 0, 'l'},
				{"write-on-image", required_argument, 0, 'm'},
				{"write-initial", required_argument, 0, 'n'},
		};

		if (print_help == 1){
			cout << "Printing help for root-segmentation-xray-2018"<< endl;

			cout << "This code has two modes: segmentation with an already preprocessed dataset, or preprocessing." << endl;

			cout << "To select preprocessing, use the flag" << endl;
			cout << "--preprocessing " << endl;
			cout << " and other arguments " << endl;
			cout << std::left << setw(42) << "--input=[STRING] "<< "Mandatory, has to be a directory." << endl;
			cout << std::left << setw(42) << "--output=[STRING] " << "Mandatory, has to be a directory." << endl;
			cout << std::left << setw(42) << "--max-threads=[INT] " << "Default is what is returned from omp_get_max_threads().  " << endl;
			cout << "All other arguments are ignored." << endl;
			cout << endl << endl;

			cout << "To select segmentation use the flag" << endl;
			cout << "--segmentation " << endl;
			cout << " and other arguments " << endl;
			cout << std::left << setw(42) << "--input=[STRING] "<< "Mandatory." << endl;
			cout << std::left << setw(42) << "--output=[STRING] " << "Mandatory." << endl;
			cout << std::left << setw(42) << "--lowerthresh=[INT] " << "lower-root-threshold; If not specified, = 0" << endl;
			cout << std::left << setw(42) << "--upperthresh=[INT] " << "upper-root-threshold; If not specified, = 255" << endl;
			cout << std::left << setw(42) << "--background=[INT] " << "background-threshold; If not specified, = 0" << endl;
			cout << std::left << setw(42) << "--nu=[FLOAT] " << "floating point value.  Default is 1.  " << endl;
			cout << std::left << setw(42) << "--grid-resolution=[INT] " << "Default is 10.  " << endl;
			cout << std::left << setw(42) << "--band-size=[INT] " << "Default is 10.  " << endl;
			cout << std::left << setw(42) << "--disregard-histogram-grid=[INT{0, 1}] " << "Default is false.  " << endl;
			cout << std::left << setw(42) << "--t=[INT] " << "Maximum value for count(x). Default is 1.  " << endl;
			cout << std::left << setw(42) << "--min-contour-size=[INT] " << "termination condition for the algo. Default is 100.  " << endl;
			cout << std::left << setw(42) << "--max-threads=[INT] " << "Default is what is returned from omp_get_max_threads().  " << endl;
			cout << std::left << setw(42) << "--write-on-image=[INT{0, 1}] " << "Option for output where segmentation is written on original slices, in color.  Default is false" << endl;
			cout << std::left << setw(42) << "--write-initial=[INT{0, 1}] " << "Option for writing the initialization masks as an image sequence." << endl;
			exit(1);
		}
		/* getopt_long stores the option index here. */
		int option_index = 0;

		opt_argument = getopt_long (argc, argv, "abcdefghijkmn",
				long_options, &option_index);

		/* Detect the end of the options. */
		if (opt_argument == -1)
			break;

		switch (opt_argument)
		{
		case 0:
			/* If this option set a flag, do nothing else now. */
			if (long_options[option_index].flag != 0)
				break;
			printf ("option %s", long_options[option_index].name);
			if (optarg)
				printf (" with arg %s", optarg);
			printf ("\n");
			break;

		case 'a':
			puts ("option -a\n");
			input_directory = optarg;
			cout << "Set " << input_directory << endl;
			break;

		case 'b':
			puts ("option -b\n");
			output_directory = optarg;
			cout << "Set " << output_directory << endl;
			break;

		case 'c':
			lower_threshold = FromString<int>(optarg);
			cout << "Lower threshold " << lower_threshold << endl;
			break;
		case 'd':
			upper_threshold = FromString<int>(optarg);
			cout << "Upper threshold " <<upper_threshold << endl;
			break;
		case 'e':
			background_threshold = FromString<int>(optarg);
			cout << "background threshold " <<background_threshold << endl;
			break;
		case 'f':
			nu = FromString<double>(optarg);
			break;
		case 'g':
			grid_resolution = FromString<int>(optarg);
			break;
		case 'h':
			band_size = FromString<int>(optarg);
			break;
		case 'i':
			disregard_histogram_grid = FromString<int>(optarg);
			break;
		case 'j':
			t = FromString<int>(optarg);
			break;
		case 'k':
			min_active_contour_size = FromString<int>(optarg);
			break;
		case 'l':
			max_number_threads = FromString<int>(optarg);
			break;
		case 'm':
			write_on_image = FromString<int>(optarg);
			break;
		case 'n':
			write_initial = FromString<int>(optarg);
			break;

		case '?':
			/* getopt_long already printed an error message. */
			break;

		default:
			abort ();
		}
	}


	if (do_preprocessing== false && do_segmentation== false){
		cout << "Options are --help, --preprocessing, or --segmentation. None were chosen.  Quitting." << endl;
		exit(1);
	}

	// then do some bad input processing.

	// are input and output dirs empty?
	// does the directory exist?
	struct stat info;

	if (output_directory.size() > 0 && input_directory.size() > 0){
		if( stat( output_directory.c_str(), &info ) != 0 ){
			cout << "cannot access " << output_directory << endl;
			exit(1);
		}
		else if( info.st_mode & S_IFDIR ){  // S_ISDIR()
			//cout << "Is a directory " << output_directory << endl;;
		}

		if( stat( input_directory.c_str(), &info ) != 0 ){
			cout << "cannot access " << input_directory << endl;
			exit(1);
		}
		else if( info.st_mode & S_IFDIR ){  // S_ISDIR()

		}
	}	else {
		if (output_directory.size() == 0){
			cout << "output directory was empty" << endl; exit(1);
		}

		if (input_directory.size() == 0){
			cout << "input directory was empty" << endl; exit(1);
		}
	}

	EnsureDirHasTrailingBackslash(input_directory);
	EnsureDirHasTrailingBackslash(output_directory);

	if (do_segmentation== true){

		// do those dirs exist?

		// s >= b?

		// bounds on thresholds?
		if (band_size <= 0){
			cout << "band size is expected to be > 0" << endl;
			exit(1);
		}
		if (grid_resolution <= 0){
			cout << "grid resolution is expected to be > 0" << endl;
			exit(1);
		}

		if (grid_resolution < band_size){
			cout << "grid resolution had to be >= band_size size.  Please consult the README." << endl;
			exit(1);
		}


		omp_set_num_threads(max_number_threads);

		string filename = output_directory + "arguments.txt";
		ofstream out;
		out.open(filename.c_str());

		out << "arguments: " << endl;
		out << "--segmentation \\" << endl;
		out << "--input=" << input_directory << " \\" << endl;
		out << "--output=" << output_directory << " \\" << endl;
		out << "--lowerthresh=" << lower_threshold << " \\" << endl;
		out << "--upperthresh=" << upper_threshold << " \\" << endl;
		out << "--background=" << background_threshold << " \\" << endl;
		out << "--nu=" << nu << " \\" << endl;
		out << "--grid-resolution=" << grid_resolution << " \\" << endl;
		out << "--band-size=" << band_size << " \\" << endl;
		out << "--disregard-histogram-grid=" << disregard_histogram_grid << " \\" << endl;
		out << "--t=" << t << " \\" << endl;
		out << "--min-contour-size=" << min_active_contour_size << " \\" << endl;
		out << "--max-threads=" << max_number_threads << " \\" << endl;
		out << "--write-on-image=" << write_on_image  << " \\" << endl;
		out << "--write-initial=" << write_initial << endl;

		out.close();


		LevelSetMainFunction(input_directory, output_directory, lower_threshold, upper_threshold, background_threshold, nu,
				grid_resolution, band_size, disregard_histogram_grid, t, min_active_contour_size, max_number_threads, write_on_image, write_initial);


	}	else {
		// do the preprocessing.
		string filename = output_directory + "arguments.txt";
		ofstream out;
		out.open(filename.c_str());
		out << "--preprocessing \\" << endl;
		out << "--input=" << input_directory << "\\" << endl;
		out << "--output=" << output_directory << "\\" << endl;
		out << "--max-threads=" << max_number_threads << "\\" << endl;
		out.close();

		vector<string> list_files;
		ReadDirectory(input_directory, list_files);
		WriteSlices(input_directory, list_files, output_directory, max_number_threads);

	}


	/* Instead of reporting ‘--verbose’
     and ‘--brief’ as they are encountered,
     we report the final status resulting from them. */
	//if (verbose_flag)
	// puts ("verbose flag is set");

	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
		printf ("non-option ARGV-elements: ");
		while (optind < argc)
			printf ("%s ", argv[optind++]);
		putchar ('\n');
	}




	exit (0);
}


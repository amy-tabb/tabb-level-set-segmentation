
#ifndef DIRECTORYFUNCTIONS_HPP_
#define DIRECTORYFUNCTIONS_HPP_

#include "Includes.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>

using std::vector;
using std::string;


void ReadDirectory(string dir_name, vector<string>& content_names);


#endif /* DIRECTORYFUNCTIONS_HPP_ */

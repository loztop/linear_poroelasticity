   // C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <sstream>

// Basic include file needed for the mesh functionality.
#include "libmesh.h"

#include "assemble.h"

using namespace std;


void read_options(unsigned int &  n_timesteps,unsigned int &  N_eles, std::string& output_file_name, std::string& result_file_name,  int& argc, char**& argv){

		n_timesteps= atoi( argv[1] );
		N_eles= atoi( argv[2] );
		output_file_name =  argv[3] ;
		result_file_name =  argv[4] ;

    std::cout<<"n_timesteps "<< n_timesteps <<" \n";
    std::cout<<"N_eles "<< N_eles <<" \n";
    std::cout<<"output_file_name "<< output_file_name <<" \n";
		std::cout<<"result_file_name "<< result_file_name <<" \n";

}

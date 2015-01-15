//
//  main.cpp
//


#include <iostream>

#include "amf.h"

////////////////////////////////////////
///////////   MAIN PROGRAM     /////////
////////////////////////////////////////

int main(int argc,char** argv)
{
	if (argc < 3)
	{
		std::cout << std::endl << "How to execute this program : \n" << std::endl
		<< "./name_executable <path_to_URM_Training> <path_to_ICM> <path_to_param_file>"
		<< "\n\n\n" ;
		return 0;
	}

	// Create an instance of the AMF solver 

	AMF amf(argv[1],argv[2],argv[3]) ;

	// Print some data

	std::cout << "We have " << amf.get_U().n_rows << " users" <<std::endl;
	std::cout << "We have " << amf.get_U().n_cols << " latent factors" <<std::endl;
	std::cout << "We have " << amf.get_V().n_rows << " labels" <<std::endl;
	std::cout << "We have " << amf.get_V().n_cols << " items" <<std::endl;
	// std::cout << "We have " << amf.get_n_max_iter_gradient() << " gradient iterations" <<std::endl;

	// Solve

	amf.solve_With_Log();







	return 0;
}

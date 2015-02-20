//
//  main_tuning.cpp
//


#include <iostream>

#include "amf.h"

////////////////////////////////////////
///////////   MAIN PROGRAM     /////////
////////////////////////////////////////

int main(int argc,char** argv)
{
	if (argc < 4)
	{
		std::cout << std::endl << "How to execute this program : \n" << std::endl
		<< "./name_executable <path_to_URM_Training> <path_to_ICM> <path_to_param_file> <path_to_urm_validation>"
		<< "\n\n\n" ;
		return 0;
	}

	// Create an instance of the AMF solver 

	AMF amf(argv[1],argv[2],argv[3]) ;

	// Solve

	amf.solve_for_tuning(argv[4]);


	// Export results

	// amf.export_Results();








	return 0;
}

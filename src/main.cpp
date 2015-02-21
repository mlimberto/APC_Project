//
//  main.cpp
//
//  Authors : Lara Bombardieri, Andrea Brandoli, Matteo Limberto


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

	// Solve

	amf.solve_With_Log();

	// Export results

	amf.export_Results();








	return 0;
}

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
		<< "./name_executable <path_to_URM> <path_to_ICM> <path_to_param_file>"
		<< "\n\n\n" ;
		return 0;
	}

	// Create an instance of the AMF solver 

	std::cout << "Creating instance of AMF ... " ;

	AMF amf(argv[1],argv[2],argv[3]) ;

	std::cout << "Done!" << std::endl;

	amf.solve();





	return 0;
}

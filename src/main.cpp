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
	std:: cout << "Hello world" << std::endl ;

	// Create instance for class AMF 

	AMF solver;

	solver.set_lambda(1.0); 

	std::cout << solver.get_lambda() << std::endl;

	Mat<double> U(5,5);

	return 0;
}

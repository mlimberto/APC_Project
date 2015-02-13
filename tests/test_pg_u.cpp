/* TEST PG_U */

#include <iostream>
#include <fstream>

#include <armadillo>

#include "my_utils.h"
#include "amf.h"

using namespace arma;
using namespace std;

bool test_PG_U(AMF & amf )
{
	// Set H,V

	// Voglio A matrice 3x5 identità

	amf.H_old_.eye(3,3);

	amf.V_old_.eye(3,5);

	(amf.H_old_ * amf.V_old_).print("A matrix");

	// Initialize U to zeros

	amf.U_.randu(5,3);
	amf.U_ = amf.U_*10.0;

	amf.U_.print("U matrix");

	// Run the test

	amf.solve_pg_U();

	amf.U_.print("U matrix");


	return true;
}

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

	std::cout << "Creating instance of AMF ... "  << std::endl;

	AMF amf(argv[1],argv[2],argv[3]) ;

	cout << ((test_PG_U(amf))? "Success " : "Failed" ) << endl;


	return 0;
}
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

	amf.H_old_.eye(3,3);
	amf.H_old_(2,1) = 1.0;
	amf.H_old_(0,2) = 1.0;

	amf.V_old_.eye(3,5);
	amf.V_old_(2,3) = 0.5;
	amf.V_old_(0,3) = 0.5;
	amf.V_old_(1,4) = 1,0;


	(amf.H_old_ * amf.V_old_).print("A matrix");

	// Set URM

	sp_mat URM(5,5);
	URM(0,0) = URM(1,1) = URM(2,2) = URM(3,3) = URM(4,4) =3.0; 
	URM(1,1) = 4.0;
	URM(3,2) = 2.0;

	URM.print("URM Matrix");

	swap(amf.URM_,URM);

	// Initialize U

	amf.U_.randu(5,3);
	amf.U_ = amf.U_*10.0;

	amf.U_.print("U matrix");

	// Initialize U_old

	amf.U_old_ = amf.U_;

	// Set toll gradient
	amf.toll_gradient_ = 1e-5;

	// Run the test

	amf.solve_pg_U();

	amf.U_.print("U matrix");

	cout << "End of the test function" << endl;

	swap(URM,amf.URM_); // Swap back to avoid memory ownership issues

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
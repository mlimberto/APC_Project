/* TEST PG_U_2 */

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

	amf.H_old_.ones(3,3);
	amf.H_old_ = amf.H_old_;



	sp_mat V_old(3,5);
	V_old(0,0) = V_old(1,1) = V_old(2,2) = 1.0;
	V_old(2,2) = V_old(1,2) = 0.5; 
	V_old(2,3) = V_old(1,4) = 1.0;

	amf.V_old_ = V_old;

	amf.H_old_.print("H matrix");
	amf.V_old_.print("V matrix");


	(amf.H_old_ * amf.V_old_).print("A matrix");

	// Set URM

	sp_mat URM(5,5);
	URM(0,0) = 0.5 ; 
	URM(0,1) = 0.5; 
	// URM(0,2) = 1.0; 
	// URM(1,0) = 1.0;
	URM(1,1) = -3.0; 
	URM(1,2) = 1.0; 
	URM(2,0) = -4; 
	// URM(2,1) = 3.0;
	URM(2,2) = 2.0; 
	URM(3,3) = URM(3,4) = 4.0;
	URM(4,4) = URM(4,2) = 1.0;


	URM.print("URM Matrix");

	amf.URM_Tr_ = URM;

	// Initialize U_old

	amf.U_old_.randu(5,3);
	amf.U_old_ = amf.U_old_*10.0;

	amf.U_old_.print("U_old matrix");

	// Set toll gradient
	amf.toll_gradient_ = 1e-5;

	// Run the test

	amf.solve_pg_U_With_Log();

	amf.U_.print("U matrix");

	(URM - amf.U_*amf.H_old_*amf.V_old_).print("S - UHV");

	cout << "Lambda = " << amf.lambda_ << endl;

	cout << "End of the test function" << endl;

	return true;
}

// bool test_PG_U(AMF & amf )
// {
// 	// Set H,V

// 	amf.H_old_.eye(3,3);
// 	amf.H_old_ = amf.H_old_*0.5;



// 	sp_mat V_old(3,3);
// 	V_old(0,0) = V_old(1,1) = V_old(2,2) = 1.0; 

// 	amf.V_old_ = V_old;

// 	amf.H_old_.print("H matrix");
// 	amf.V_old_.print("V matrix");


// 	(amf.H_old_ * amf.V_old_).print("A matrix");

// 	// Set URM

// 	sp_mat URM(3,3);
// 	URM(0,0) = 0.5 ; 
// 	URM(0,1) = 0.5; 
// 	URM(0,2) = 1.0; 
// 	URM(1,0) = 1.0;
// 	URM(1,1) = -3.0; 
// 	URM(1,2) = 1.0; 
// 	URM(2,0) = -4; 
// 	URM(2,1) = 3.0;
// 	URM(2,2) = 2.0; 


// 	URM.print("URM Matrix");

// 	amf.URM_Tr_ = URM;

// 	// Initialize U_old

// 	amf.U_old_.randu(3,3);
// 	amf.U_old_ = amf.U_old_*10.0;

// 	amf.U_old_.print("U_old matrix");

// 	// Set toll gradient
// 	amf.toll_gradient_ = 1e-5;

// 	// Run the test

// 	amf.solve_pg_U_With_Log();

// 	amf.U_.print("U matrix");

// 	(URM - amf.U_*amf.H_old_*amf.V_old_).print("S - UHV");

// 	cout << "Lambda = " << amf.lambda_ << endl;

// 	cout << "End of the test function" << endl;

// 	return true;
// }


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
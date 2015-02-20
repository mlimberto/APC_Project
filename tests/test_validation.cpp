/* TEST PG_U_2 */

#include <iostream>
#include <fstream>

#include <armadillo>

#include "my_utils.h"
#include "amf.h"

using namespace arma;
using namespace std;

bool test_validation(AMF & amf)
{
	// Create a URM_Validation matrix

	umat URM_Locations; 
	URM_Locations << 0 << 1 << 2 << 0 << endr 
				  << 0 << 1 << 2 << 1 << endr;

	vec values;
	values << 2 << 2 << 2 << 5 << endr;

	sp_mat URM_Val(URM_Locations,values);

	URM_Val.print("URM validation");

	amf.URM_Val_Location_Matrix_ = URM_Locations ;
	amf.URM_Val_Values_ = values;

	// Set U H V 

	mat U(3,3,fill::ones);
	mat H(3,3,fill::eye);
	sp_mat V(3,3);
	V(0,0) = V(1,1) = V(2,2) = 1.0;

	amf.U_old_ = U;
	amf.H_old_ = H;
	amf.V_old_ = V;

	std::cout << amf.evaluate_Against_URM_Validation() << std::endl;

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

	cout << ((test_validation(amf))? "Success " : "Failed" ) << endl;


	return 0;
}
#include <iostream>

#include "amf.h"

using namespace std;

int main(int argc,char** argv)
{
	cout << "Testing initialization of parameters" << endl;

	AMF amf;

	amf.initialize_Parameters("test_parameters.txt") ;

	cout << "Lambda has been set to " << amf.get_lambda() << endl;

	cout << "Lat factors has been set to " << amf.get_n_latent_factors() << endl;

	cout << "Toll has been set to " << amf.get_toll() << endl;

	cout << "n_max_iter has been set to " << amf.get_n_max_iter() << endl;

	cout << "n_max_iter_gradient has been set to " << amf.get_n_max_iter_gradient() << endl;

	cout << "toll_gradient has been set to " << amf.get_toll_gradient() << endl;

	cout << "gradient_step has been set to " << amf.get_gradient_step() << endl;

	cout << ( (amf.get_lambda()==1.5) ? "Succeed" : "Failed") << endl;
	

	return 0;
}
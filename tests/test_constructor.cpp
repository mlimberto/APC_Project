#include <iostream>

#include "amf.h"

using namespace std;

int main(int argc,char** argv)
{

	cout << argc << endl;

	for (auto i = 0 ; i < argc ; ++i)
		cout << argv[i] << endl;



	cout << "Testing constructor of an instance of AMF" << endl;

	AMF amf(argv[1],argv[2],argv[3]);


	cout << amf.get_lambda() << endl;
	cout << amf.get_n_max_iter() << endl;
	

	return 0;
}
#include <iostream>
#include <armadillo>

#include "my_utils.h"

using namespace std;
using namespace arma;

int main(int argc,char **argv)
{
	SpMat<double> V(3,3);
	SpMat<int> X(3,3);

	V(0,0) = 0.5;
	X(0,0) = 1;
	V(1,0) = 0.5;
	// X(1,0) = 1;

	V(2,2) = 1.0;
	X(2,2) = 1;

	V(1,1) = 1.0;
	X(1,1) = 1;

	X(0,2) = 1;

	cout << V << X;

	bool ans = check_V_Constraint(V,X) ;

	cout <<  ( ans? "Success" : "Failed" ) << endl;

	X(1,0) = 1;

	cout <<  ( check_V_Constraint(V,X)? "Success" : "Failed" ) << endl;

	return 0;
}
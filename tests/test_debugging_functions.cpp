#include <iostream>
#include <armadillo>

#include "my_utils.h"

using namespace std;
using namespace arma;

int main(int argc,char **argv)
{

	// Test check_V_Constraint

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

	// Test check_Positive_Matrix
	cout << "Testing positivity of matrix" << endl;

	mat U = randu<mat>(4,5);
	cout << U;

	cout <<  ( check_Positive_Matrix(U)? "Success" : "Failed" ) << endl;	

	U(0,0) = -0.1;
	cout << U;
	
	cout <<  ( check_Positive_Matrix(U)? "Success" : "Failed" ) << endl;

	return 0;
}
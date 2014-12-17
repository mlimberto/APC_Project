#include <iostream>

#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc,char **argv)
{
	mat A = randu<mat>(4,5);
	mat B = randu<mat>(4,5);

	// cout << A*B.t() << endl;
	cout << A << endl;

	B = cumsum(A);

	cout << B << endl;

	SpMat<double> V(5,5);

	cout << V << endl; 

	return 0;
}
#include <iostream>

#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc,char **argv)
{

	// Test iterators full matrix

	mat A = randu<mat>(4,5);
	A.ones(); // Set all the elements of A to 1 

	// Using AUTO
	double sum = 0;

	for (auto i = A.begin_col(1); i != A.end_col(1) ; ++i)
	{
		sum += *i ;
	}

	cout << sum << endl;

	// Using COL_ITERATOR
	sum = 0 ;

	for (mat::col_iterator i = A.begin_col(1) ; i != A.end_col(1) ; ++i)
	{
		sum += *i ;
	}

	cout << sum << endl;

	// Using CONST_COL_ITERATOR to enforce read-only access
	sum = 0 ;

	for (mat::const_col_iterator i = A.begin_col(1) ; i != A.end_col(1) ; ++i)
	{
		sum += *i ;
	}

	cout << sum << endl;


	// ITERATOR ON SPARSE MATRIX

	SpMat<double> V(5,5);

	V(0,0) = 4.0;
	V(4,1) = 3.0;
	V(2,2) = 1.0;

	cout << V ; 

	// Using AUTO
	sum = 0 ;

	for (auto i = V.begin_col(1) ; i != V.end_col(1) ; ++i)
	{
		sum += *i ;
	}

	cout << sum << endl;

	// Using library-defined type

	sum = 0 ;

	for (sp_mat::const_iterator i = V.begin_col(1) ; i != V.end_col(1) ; ++i)
	{
		sum += *i ;
	}

	cout << sum << endl;

	// Does the same work using an integer matrix? 

	SpMat<int> X(5,5);

	X(0,0) = 4;
	X(4,1) = 3;
	X(2,2) = 1;

	cout << X ; 

	// Using AUTO
	int int_sum = 0 ;

	for (auto i = X.begin_col(1) ; i != X.end_col(1) ; ++i)
	{
		int_sum += *i ;
	}

	cout << int_sum << endl;

	// Using library-defined type

	int_sum = 0 ;

	for (SpMat<int>::const_iterator i = X.begin_col(1) ; i != X.end_col(1) ; ++i)
	{
		int_sum += *i ;
	}

	cout << int_sum << endl;

	return 0;
}
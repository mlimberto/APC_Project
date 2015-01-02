//
//  my_utils.cpp
//  
//
//
//

#include "my_utils.h"

using std::string;
using std::vector;
using std::ifstream;




double build_S(int i, int j, SpMat<double> URM, Mat<double> U, Mat<double> H, SpMat<double> V )
{

	return 0;
}


void get_Positive_Matrix(Mat<double> &U)
{
	for (auto i = U.begin() ; i != U.end() ; ++i)
	{
		if (*i < 0)
			*i = 0;
	}
}

///////////////////////////////////
////// DEBUGGING FUNCTIONS ////////
///////////////////////////////////


bool check_V_Constraint(const SpMat<double>& V,const SpMat<int>& X)
{
	// CHECK 1
	// All zero elements of X must be zero elements of V 

	if (V.n_cols != X.n_cols || V.n_rows != X.n_rows)
	{
		std::cerr << "Inconsistent matrices !!!" << std::endl;
	}

	// I must loop over all the zero elements of X (it's not enough to loop
	// over the non-zero elements of X unfortunately)
	// TODO This loop is inefficient, we could use the CSC format to perform
	// this check, not very important though

	for (auto i = 0 ; i<X.n_cols ; ++i)
	{
		for (auto j = 0 ; j < X.n_rows ; ++j) 
		{	
			if (X(i,j)==0 && V(i,j)!=0)
			{
				#ifndef NDEBUG
				std::cerr << "Mismatch found at position (" << i << "," << j << ") " << std::endl;
				#endif
				return false;
			}
		}
	}

	#ifndef NDEBUG
	std::cout << "Consistency with ICM passed!" << std::endl;
	#endif

	// CHECK 2
	// For every column the sum of the elements must be equal to 1
	double eps = 1e-4 ;

	// Loop on each column
	for (std::size_t n=0 ; n< V.n_cols ; n++)
	{
		double sum(0);

		// Loop on the iterator for element access

		for (SpMat<double>::const_iterator i = V.begin_col(n) ; i != V.end_col(n) ; ++i)
		{
			sum += *i ;
		}


		if (sum > 1 + eps || sum < 1 - eps)
		{
			#ifndef NDEBUG
			std::cerr << "Sum is not consistent on column " << n << std::endl;
			#endif
			return false;
		}

	}

	#ifndef NDEBUG
	std::cout << "Consistency with probabilistic interpretation passed! " << std::endl;
	#endif

	return true;
}


bool check_Positive_Matrix(const Mat<double> &U)
{
	// Loop on the elements of the matrix to check positivity 

	for (Mat<double>::const_iterator i = U.begin() ; i != U.end() ; ++i )
	{
		if ( (*i)<0 )
		{
			#ifndef NDEBUG
			std::cerr << "Found negative element!" << std::endl;
			#endif
			return false;
		}
	}

	return true;
}



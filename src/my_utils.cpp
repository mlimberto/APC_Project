//
//  my_utils.cpp
//  
//
//
//

#include "my_utils.h"

using std::string;
using std::vector;


// ricordarsi che gli indici di riga e colonna in c++ partono da 0 !!

double build_S(uword i, uword j,const SpMat<double>& URM,const Mat<double>& U,const Mat<double>& H,const SpMat<double>& V )
{
	if (URM(i,j) != 0){
		return URM(i,j);
	}
	else
	{
		return as_scalar(U.row(i)*H*V.col(j));
	}
}

double evaluate_Obj_Function(const SpMat<double>& URM,const Mat<double>& U,
							 const Mat<double>& H,const SpMat<double>& V,
							 const Mat<double>& U_old,const Mat<double>& H_old,
							 const SpMat<double> V_old, const double lambda)
{

	double obj(0);

	for ( size_t i(0), i < URM.n_rows , i++){
		for ( size_t j(0), j < URM.n_cols, j++){
			obj += pow( build_S( i, j, URM, U_old, H_old, V_old) - U.row(i)*H*V.col(j) , 2 );
		}
	}

	return as_scalar( obj + lambda * ( pow( norm (U,"fro"), 2) + pow( norm (H,"fro"), 2) ) );
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



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

bool read_URM_From_File(vector<unsigned int> &rows , 
						vector<unsigned int> &cols , 
						vector<double> &values ,
						string matrix_filename ) 
{
	// Check for the vectors to be empty and if not clear them
	if (!rows.empty()) rows.clear();
	if (!cols.empty()) cols.clear();
	if (!values.empty()) values.clear();

	ifstream matrix_file(matrix_filename);
	unsigned int ri,ci;
	double val;

	while (matrix_file >> ri >> ci >> val)
	{
		rows.push_back(ri);
		cols.push_back(ci);
		values.push_back(val);
	}

	return true;
}

bool read_ICM_From_File(vector<unsigned int> &rows , 
						vector<unsigned int> &cols , 
						vector<unsigned int> &values ,
						string matrix_filename  ) 
{
	// Check for the vectors to be empty and if not clear them
	if (!rows.empty()) rows.clear();
	if (!cols.empty()) cols.clear();
	if (!values.empty()) values.clear();

	ifstream matrix_file(matrix_filename);
	unsigned int ri,ci,val;

	while (matrix_file >> ri >> ci >> val)
	{
		rows.push_back(ri);
		cols.push_back(ci);
		values.push_back(val);
	}

	return true;
}


double build_S(int i, int j, SpMat<double> URM, Mat<double> U, Mat<double> H, SpMat<double> V )
{

	return 0;
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



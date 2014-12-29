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

	return false;
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


bool check_V_Constraint(SpMat<double> V)
{

	return true;
}



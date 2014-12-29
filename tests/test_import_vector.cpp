#include <iostream> 
#include <string>
#include <vector>

#include "my_utils.h"

int main(int argc,char** argv)
{

	std::vector<unsigned int> rows,cols,values_ICM;
	std::vector<double> values;

	std::string matrix_filename("../dataset/Large_Dataset/urm_converted.txt");

	std::cout << "Testing input/output routines. " << std::endl;

	// TEST URM

	std::cout << "Testing URM ..." << std::endl;

	read_URM_From_File(rows,cols,values,matrix_filename);

	std::cout <<( (rows.size() == 8862158) ? "Test for rows passed" : "Test for rows failed" )<< std::endl;
	std::cout <<( (cols.size() == 8862158) ? "Test for cols passed" : "Test for cols failed" )<< std::endl;
	std::cout <<( (values.size() == 8862158) ? "Test for values passed" : "Test for values failed" )<< std::endl;


	// TEST ICM

	std::cout << "Testing ICM ..." << std::endl;

	matrix_filename = "../dataset/Large_Dataset/icm_converted.txt" ;

	read_ICM_From_File(rows,cols,values_ICM,matrix_filename);

	std::cout <<( (rows.size() == 383352) ? "Test for rows passed" : "Test for rows failed" )<< std::endl;
	std::cout <<( (cols.size() == 383352) ? "Test for cols passed" : "Test for cols failed" )<< std::endl;
	std::cout <<( (values_ICM.size() == 383352) ? "Test for values passed" : "Test for values failed" )<< std::endl;


	return 0;
}
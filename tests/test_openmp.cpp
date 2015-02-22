//
//  main.cpp
//
//  Authors : Lara Bombardieri, Andrea Brandoli, Matteo Limberto


#include <iostream>
#include <chrono>

#include "amf.h"

////////////////////////////////////////
///////////   MAIN PROGRAM     /////////
////////////////////////////////////////

bool test_openmp(AMF &amf)
{
	#ifdef AMFOPENMP
	std::cout<< "Using OPENMP ... " << std::endl;
	#endif

	// Evaluate objective function ... 

	double res(0);

	auto begin = std::chrono::high_resolution_clock::now();

	res = evaluate_Obj_Function(amf.URM_Tr_ , amf.U_old_ ,amf.H_old_,amf.V_old_ ,
						amf.U_old_ ,amf.H_old_, amf.V_old_ , 1.0 ) ;

	auto end = std::chrono::high_resolution_clock::now();
	
	std::cout <<"Objective function : " << res << std::endl;
	std::cout <<"Evaluation took "<< std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << "ms" << std::endl;
	
	// Evaluate linear part of gradient for V ...

	// amf.U_ = amf.U_old_;
	// amf.H_ = amf.H_old_;
	amf.n_max_iter_gradient_ = 1;
	// amf.solve_V();


	// Evaluate linear part of gradient for H ...
	// amf.solve_pg_H();
	amf.solve_pg_U_With_Log();
	amf.solve_pg_U();


	return true;
}

int main(int argc,char** argv)
{
	if (argc < 3)
	{
		std::cout << std::endl << "How to execute this program : \n" << std::endl
		<< "./name_executable <path_to_URM_Training> <path_to_ICM> <path_to_param_file>"
		<< "\n\n\n" ;
		return 0;
	}

	// Create an instance of the AMF solver 

	AMF amf(argv[1],argv[2],argv[3]) ;

	test_openmp(amf);

	return 0;
}

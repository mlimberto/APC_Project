//
//  proj_gradient.cpp
//  
//
//
//

#include "amf.h"


//////////////////////////
////  METHODS       //////
//////////////////////////


void AMF::solve_pg_U()
{
	// Allocate and compute the temporary matrix A := H*V
	// and the gradient matrix G
	#ifndef NDEBUG
	std::cout << "Solving projected gradient for U " << std::endl;
	#endif

	mat G(U_.n_rows,U_.n_cols,fill::zeros);

	mat A = H_old_*V_old_;


	// Run the loop

	bool stop_criterion = false;

	for (unsigned int n=0 ; (n < n_max_iter_gradient_ ) && (!stop_criterion) ; ++n ) 
	{
		// Perform the gradient step 
		solve_pg_U_One_Iteration(G,A);

		// Evaluate stop criterion

	}


}

void AMF::solve_pg_U_With_Log()
{

}

void AMF::solve_pg_U_One_Iteration(mat &G, mat &A)
{
	// Compute the gradient


	// Update U_ 
	U_ = U_ - gradient_step_*G ;

	// Projection step
	get_Positive_Matrix(U_);


}

void AMF::solve_pg_H()
{

}

void AMF::solve_pg_H_With_Log()
{

}

void AMF::solve_pg_H_One_Iteration(mat &G,mat &A)
{

}
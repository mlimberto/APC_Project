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
	#ifndef NDEBUG
	std::cout << "Solving projected gradient for U " << std::endl;
	#endif

	// Allocate and compute the temporary matrix A := H*V
	// and the gradient matrix G

	mat G(U_.n_rows,U_.n_cols,fill::zeros);

	mat A = H_old_*V_old_; // This matrix shall be computed only once and stored

	// Run the loop

	bool stop_criterion = false;

	for (unsigned int n=0 ; (n < n_max_iter_gradient_ ) && (!stop_criterion) ; ++n ) 
	{
		// #ifndef NDEBUG
		// std::cout <<"Iteration " << n << std::endl;
		// #endif


		// Perform the gradient step 
		solve_pg_U_One_Iteration(G,A);

		// Evaluate stop criterion


	}


}

void AMF::solve_pg_U_With_Log()
{

}

void AMF::solve_pg_U_One_Iteration(mat &G,const mat &A)
{
	// Compute the gradient

	for (uword x =0 ; x< U_.n_rows ; ++x)
	{
		for (uword y = 0 ; y < U_.n_cols ; ++y)
		{
			double qq = 0;

			// Parte quadratica
			for (uword k = 0 ; k < U_.n_cols ; ++k)
				for(uword j = 0 ; j < U_.n_cols ; ++j)
					qq += U_(x,k)*A(k,j)*A(y,j);
				
			// Parte lineare
			double ll = 0 ;

			for (uword k = 0 ; k < A.n_cols ; ++k)
				ll += build_S(x,k,URM_,U_old_,H_old_, V_old_) * A(y,k);

			// Aggiorna la matrice

			G(x,y) = 2*qq - 2*ll ; // La parte legata all'overfitting viene aggiunta dopo

			}
		}


	// Update U_ 
	U_ = U_ - gradient_step_*G -2*gradient_step_*lambda_*U_ ;

	// #ifndef NDEBUG
	// std::cout << std::scientific << norm(U_,"fro")  << std::endl;
	// #endif

	// Projection step
	get_Positive_Matrix(U_);

	// #ifndef NDEBUG
	// std::cout << std::scientific << norm(U_,"fro") << std::endl;
	// #endif

	// #ifndef NDEBUG
	// U_.print("U matrix");
	// G.print("G matrix");
	// #endif	


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
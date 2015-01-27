//
//  proj_gradient.cpp
//  
//
//
//

#include "amf.h"

using namespace arma;

//////////////////////////
////  METHODS       //////
//////////////////////////


void AMF::solve_pg_U()
{
	#ifndef NDEBUG
	std::cout << "Solving projected gradient for U " << std::endl;
	#endif
	
	// Initialize U (check what is best)
	U_ = U_old_;

	// Allocate and compute the temporary matrix A := H*V
	// and the gradient matrix G
	mat G(U_.n_rows,U_.n_cols,fill::zeros);
	mat A = H_old_*V_old_; // This matrix shall be computed only once and stored

	// Run the loop

	bool stop_criterion = false;

	double prec_obj(0) , curr_obj(0);

	for (unsigned int n=0 ; (n < n_max_iter_gradient_ ) && (!stop_criterion) ; ++n ) 
	{
		// #ifndef NDEBUG
		// std::cout <<"Gradient method for U, iteration " << n << std::endl;
		// #endif


		// Perform the gradient step 
		solve_pg_U_One_Iteration(G,A);

		// Evaluate stop criterion
		prec_obj = curr_obj;
		curr_obj = evaluate_Obj_Function(URM_Tr_,U_,H_old_,V_old_,U_old_,H_old_,V_old_,lambda_);

		if (abs(curr_obj - prec_obj)/curr_obj < toll_gradient_)
			stop_criterion = true;

	}


}

void AMF::solve_pg_U_With_Log()
{
	std::ofstream logfile;
	logfile.open("log_pg_u.txt");

	// Initialize U (check what is best)
	U_ = U_old_;

	mat G(U_.n_rows,U_.n_cols,fill::zeros);
	mat A = H_old_*V_old_;

	bool stop_criterion = false;

	double prec_obj(0) , curr_obj(0);

	for (unsigned int n=0 ; (n < n_max_iter_gradient_ ) && (!stop_criterion) ; ++n ) 
	{
		solve_pg_U_One_Iteration(G,A);

		// Evaluate stop criterion
		prec_obj = curr_obj;
		curr_obj = evaluate_Obj_Function(URM_Tr_,U_,H_old_,V_old_,U_old_,H_old_,V_old_,lambda_);

		// if (abs(curr_obj - prec_obj)/curr_obj < toll_gradient_)
		// 	stop_criterion = true;

		// Save information on logfile
		logfile << curr_obj << "\n";

		// Print information
		std::cout << "Iteration " << n+1 << " : Objective function = " << curr_obj << std::endl;

	}

	#ifndef NDEBUG
	U_.print("U matrix");
	#endif

	logfile.close();
}

void AMF::solve_pg_U_One_Iteration(mat &G,const mat &A)
{
	// Compute the gradient
	for (uword x =0 ; x< U_.n_rows ; ++x)
	{
		for (uword y = 0 ; y < U_.n_cols ; ++y)
		{
			// Parte quadratica
			double qq = 0;
			for (uword k = 0 ; k < U_.n_cols ; ++k)
				for(uword j = 0 ; j < U_.n_cols ; ++j)
					qq += U_(x,k)*A(k,j)*A(y,j);
				
			// Parte lineare
			double ll = 0 ;
			for (uword k = 0 ; k < A.n_cols ; ++k)
				ll += build_S(x,k,URM_Tr_,U_old_,H_old_, V_old_) * A(y,k);

			// Aggiorna la matrice
			G(x,y) = 2*qq - 2*ll ; // La parte legata all'overfitting viene aggiunta dopo

			}
		}




	// Update U_ 
	U_ = U_ - gradient_step_*G -2*gradient_step_*lambda_*U_ ;

	// Projection step
	get_Positive_Matrix(U_);

}

void AMF::solve_pg_H(){
    #ifndef NDEBUG
    std::cout << "Solving projected gradient for H " << std::endl;
    #endif

    // Initialize U (check what is best)
    H_ = H_old_;

    // Allocate and compute the temporary matrix A := H*V
    // and the gradient matrix G
    mat G(H_.n_rows,H_.n_cols,fill::zeros);

    // Run the loop

    bool stop_criterion = false;

    double prec_obj(0) , curr_obj(0);

    for (unsigned int n=0 ; (n < n_max_iter_gradient_ ) && (!stop_criterion) ; ++n )
    {
        // #ifndef NDEBUG
        // std::cout <<"Gradient method for H, iteration " << n << std::endl;
        // #endif


        // Perform the gradient step
        solve_pg_H_One_Iteration(G);

        // Evaluate stop criterion
        prec_obj = curr_obj;
        curr_obj = evaluate_Obj_Function(URM_Tr_,U_,H_,V_old_,U_old_,H_old_,V_old_,lambda_);

        if (abs(curr_obj - prec_obj)/curr_obj < toll_gradient_)
            stop_criterion = true;

    }
}

void AMF::solve_pg_H_With_Log(){

    std::ofstream logfile;
    logfile.open("log_pg_h.txt");

    // Initialize H (check what is best)
    H_ = H_old_;

    mat G(U_.n_rows,U_.n_cols,fill::zeros);
    mat A = H_old_*V_old_;

    bool stop_criterion = false;

    double prec_obj(0) , curr_obj(0);

    for (unsigned int n=0 ; (n < n_max_iter_gradient_ ) && (!stop_criterion) ; ++n )
    {
        solve_pg_H_One_Iteration(G);

        // Evaluate stop criterion
        prec_obj = curr_obj;
        curr_obj = evaluate_Obj_Function(URM_Tr_,U_,H_,V_old_,U_old_,H_old_,V_old_,lambda_);

        // if (abs(curr_obj - prec_obj)/curr_obj < toll_gradient_)
        // 	stop_criterion = true;

        // Save information on logfile
        logfile << curr_obj << "\n";

        // Print information
        std::cout << "Iteration " << n+1 << " : Objective function = " << curr_obj << std::endl;

    }

    #ifndef NDEBUG
    H_.print("H matrix");
    #endif

    logfile.close();
}

void AMF::solve_pg_H_One_Iteration(mat &G){

    // Compute the gradient
    for (uword alpha =0 ; alpha< H_.n_rows ; ++alpha){
        for (uword beta = 0 ; beta < H_.n_cols ; ++beta){

            vec q(U_.n_rows,fill::zeros);

            for(uword i=0; i<q.n_elem; ++i)
                for(uword j=0; j<V_.n_cols; ++j)
                    q(i)=q(i)+(U_.row(i)*H_*V_old_.col(j)-build_S(i,j,URM_Tr_,U_old_,H_old_, V_old_))*V_old_(beta,j);

            G(alpha,beta)=2*dot(U_.col(alpha),q)+2*lambda_*H_;

            }
        }


    // Update H_
    H_ = H_ - gradient_step_*G;

    // Projection step
    get_Positive_Matrix(H_);
}

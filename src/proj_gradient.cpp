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


void AMF::solve_pg_U() // STILL WORK IN PROGRESS!!!!
{
	#ifndef NDEBUG
	std::cout << "Solving projected gradient for U " << std::endl;
	#endif
	
	// Initialize U (check what is best)
	U_ = U_old_;

	// Allocate and compute the temporary matrix A := H*V

	mat A = H_old_*V_old_; 

	// and the gradient matrix G
	mat G(U_.n_rows,U_.n_cols,fill::zeros);



	// Run the loop

	bool stop_criterion = false;

	double prec_obj(0) , curr_obj(0);

	for (unsigned int n=0 ; (n < n_max_iter_gradient_ ) && (!stop_criterion) ; ++n ) 
	{
		// #ifndef NDEBUG
		// std::cout <<"Gradient method for U, iteration " << n << std::endl;
		// #endif


		// Perform the gradient step 
		// solve_pg_U_One_Iteration(G,A);

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

	// Initialize U (with "warm-up")
	U_ = U_old_;

	mat A = H_old_*V_old_;
	mat AAt = A*(A.t() );

	mat G(U_.n_rows,U_.n_cols,fill::zeros);

	gradient_step_ = 1.0;

	// Compute the linear part of the gradient
	for (uword x =0 ; x< U_.n_rows ; ++x)
		for (uword y = 0 ; y < U_.n_cols ; ++y)
		{				
			double ll = 0 ;

			for (uword k = 0 ; k < A.n_cols ; ++k)
				ll += build_S(x,k,URM_Tr_,U_old_,H_old_, V_old_) * A(y,k);

			G(x,y) = - ll ; 
		}


	bool stop_criterion = false;

	double prec_obj(0) , curr_obj(0);

	for (unsigned int n=0 ; (n < n_max_iter_gradient_ ) && (!stop_criterion) ; ++n ) 
	{
		solve_pg_U_One_Iteration(G,A,AAt);

		// Evaluate stop criterion
		prec_obj = curr_obj;
		curr_obj = evaluate_Obj_Function(URM_Tr_,U_,H_old_,V_old_,U_old_,H_old_,V_old_,lambda_);

		if (abs(curr_obj - prec_obj)/curr_obj < toll_gradient_)
			stop_criterion = true;

		#ifndef NDEBUG
		std::cout << "old_obj " << prec_obj << " new_obj" << curr_obj << std::endl;
		#endif

		// Save information on logfile
		logfile << curr_obj << "\n";
		total_logfile_ << curr_obj << "\n";

		// Print information
		std::cout << "Iteration " << n+1 << " : Objective function = " << curr_obj << std::endl;

	}


	logfile.close();
}

void AMF::solve_pg_U_One_Iteration(mat G,const mat &A, const arma::mat &AAt)
{
	// Update quadratic part of the gradient
	G += U_*AAt;

	// Find a feasible step
	double sigma = 0.01;
	double beta = 0.1;
	bool is_feasible = false;

	// Make a first try ...
	mat U_cand = U_ - gradient_step_*G -gradient_step_*lambda_*U_ ;
	get_Positive_Matrix(U_cand);
	mat D = U_cand - U_ ;
	
	double res = (1-sigma)*dot(G,D) + 
				0.5*dot(D , D *(AAt + lambda_*eye(G.n_cols,G.n_cols) )) ;

	std::cout << "Step = " << gradient_step_ << " Value = " << res << " ";
	std::cout <<( (res <= 0)?("First is feasible so I increase the step"):("First is not feasible") )<< std::endl;

	if (res <=0)
	{
		gradient_step_ = gradient_step_ / beta;
	}
	else
		gradient_step_ = gradient_step_ * beta;

	// Let's keep on trying ...

	while(!is_feasible)
	{
		U_cand = U_ - gradient_step_*G -gradient_step_*lambda_*U_ ;
		get_Positive_Matrix(U_cand);
		D = U_cand - U_ ;

		res = (1-sigma)*dot(G,D) + 
				0.5*dot(D , D *(AAt + lambda_*eye(G.n_cols,G.n_cols) )) ;

		std::cout << "Step = " << gradient_step_ << " Value = " << res << " ";
		std::cout <<( (res <= 0)?("Feasible"):("Not feasible") )<< std::endl;

		if (res <=0)
			is_feasible = true;
		else 
			gradient_step_ = beta*gradient_step_;
	}

	std::cout << "Selected step is " << gradient_step_ << std::endl;

	// Update U_ 
	U_ = U_ - gradient_step_*G -gradient_step_*lambda_*U_ ;
	get_Positive_Matrix(U_);


}

void AMF::solve_pg_H(){
    #ifndef NDEBUG
    std::cout << "Solving projected gradient for H " << std::endl;
    #endif

    // Initialize H (with warm up)
    H_ = H_old_;

    // Allocate and compute the temporary matrix A := H*V
    // and the gradient matrix G
    mat G(H_.n_rows,H_.n_cols,fill::zeros);

    // Run the loop

    bool stop_criterion = false;

    double prec_obj(0) , curr_obj(0);

    for (unsigned int n=0 ; (n < n_max_iter_gradient_ ) && (!stop_criterion) ; ++n )
    {



        // Perform the gradient step
        // solve_pg_H_One_Iteration(G);

        // Evaluate stop criterion
        prec_obj = curr_obj;
        curr_obj = evaluate_Obj_Function(URM_Tr_,U_,H_,V_old_,U_old_,H_old_,V_old_,lambda_);

        // Added extra security measure to stop criterium
        if (abs(prec_obj - curr_obj)/curr_obj < toll_gradient_)
            stop_criterion = true;

    }
}

void AMF::solve_pg_H_With_Log(){

    std::ofstream logfile("log_pg_h.txt");

    // Initialize H (with warm up)
    H_ = H_old_;

    // Precompute UtU and VVt
	std::cout << "Compiting UtU ..." << std::endl;

    mat UtU = (U_.t())*U_;


    std::cout << "Computing VVt ..." << std::endl;

    mat VVt(V_old_.n_rows,V_old_.n_rows,fill::zeros);
	for (uword i = 0 ; i < V_old_.n_rows ; ++i )
		for (uword j = 0 ; j < V_old_.n_rows ; ++j )
		{
			for (uword k = 0 ; k < V_old_.n_cols ; ++k )
				VVt(i,j) += V_old_(i,k)*V_old_(j,k);
		}

	// Initialize gradient step
	gradient_step_ = 1.0;

    // Precompute linear part of the gradient
    mat G(H_.n_rows,H_.n_cols,fill::zeros);

    std::cout << "Computing linear part of gradient ..." << std::endl;

	for (uword k = 0 ; k < U_.n_rows ; ++k)
 	for (uword l = 0 ; l < V_old_.n_cols ; ++l)
 	{
 		double aux = build_S(k,l,URM_Tr_,U_old_,H_old_, V_old_);
		for (uword x =0 ; x< H_.n_rows ; ++x)
		for (uword y = 0 ; y < H_.n_cols ; ++y)
 				G(x,y) = G(x,y) - U_(k,x)*aux*V_old_(y,l) ;
       	
    }
    
    bool stop_criterion = false;

    double prec_obj(0) , curr_obj(0);

    for (unsigned int n=0 ; (n < n_max_iter_gradient_ ) && (!stop_criterion) ; ++n )
    {

        solve_pg_H_One_Iteration(G,UtU,VVt);

        // Evaluate stop criterion
        prec_obj = curr_obj;
        curr_obj = evaluate_Obj_Function(URM_Tr_,U_,H_,V_old_,U_old_,H_old_,V_old_,lambda_);

        if (abs(prec_obj - curr_obj)/curr_obj < toll_gradient_)
        	stop_criterion = true; 

        // Save information on logfile
        logfile << curr_obj << std::endl;
		total_logfile_ << curr_obj << "\n";


        // Print information
        std::cout << "Iteration " << n+1 << " : Objective function = " << curr_obj << std::endl;
        std::cout << abs(curr_obj - prec_obj)/curr_obj << std::endl;

    }


    logfile.close();
}

void AMF::solve_pg_H_One_Iteration(mat G,const mat &UtU, const mat &VVt){

	// Update gradient with non-constant part

    G = G + UtU*H_*VVt + lambda_*H_;


	// Find a feasible step
	double sigma = 0.01;
	double beta = 0.1;
	bool is_feasible = false;

	// Make a first try...

	mat H_cand = H_ - gradient_step_*G ;
	get_Positive_Matrix(H_cand);
	mat D = H_cand - H_ ;

	double hessian_part = dot(D,UtU*D*(VVt.t())) + lambda_*dot(D,D) ;
	hessian_part = 0.5*hessian_part;

	double res = (1-sigma)*dot(G,D) + hessian_part ;

	std::cout << "Step = " << gradient_step_ << " Value = " << res << " ";
	std::cout <<( (res <= 0)?("First is feasible so I increase the step"):("Fist is not feasible") )<< std::endl;

	if (res <=0)
		gradient_step_ = gradient_step_/beta;	
	else 
		gradient_step_ = beta*gradient_step_;

	// Let's keep on trying ...

	while(!is_feasible)
	{
		mat H_cand = H_ - gradient_step_*G ;
		get_Positive_Matrix(H_cand);
		mat D = H_cand - H_ ;

		double hessian_part = dot(D,UtU*D*(VVt.t())) + lambda_*dot(D,D) ;
		hessian_part = 0.5*hessian_part;

		double res = (1-sigma)*dot(G,D) + hessian_part ;

		std::cout << "Step = " << gradient_step_ << " Value = " << res << " ";
		std::cout <<( (res <= 0)?("Feasible"):("Not feasible") )<< std::endl;

		if (res <=0)
			is_feasible = true;
		else 
			gradient_step_ = beta*gradient_step_;
	}

	std::cout << "Selected step is " << gradient_step_ << std::endl;

    // Update H_
    H_ = H_ - gradient_step_*G;
    get_Positive_Matrix(H_);

}

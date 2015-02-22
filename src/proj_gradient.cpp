//
//  proj_gradient.cpp
//  
//
//
//  Authors : Lara Bombardieri, Andrea Brandoli, Matteo Limberto

#include "amf.h"

using namespace arma;

//////////////////////////
////  METHODS       //////
//////////////////////////


void AMF::solve_pg_U() 
{
	U_ = U_old_;

	mat A = H_old_*V_old_;
	mat AAt = A*(A.t() );

	mat G(U_.n_rows,U_.n_cols,fill::zeros);

	gradient_step_ = 1.0;

	// Compute the constant part of the gradient
	#ifndef NDEBUG
	std::cout << "Computing linear part of gradient ... " ;
	#endif

	#ifdef AMFTIME 
	auto begin = std::chrono::high_resolution_clock::now();
	#endif

	#ifdef AMFOPENMP
	#pragma omp parallel for
	#endif
	for (uword k = 0 ; k < A.n_cols ; ++k)
		for (uword x =0 ; x< U_.n_rows ; ++x)
		{
			double aux = build_S(x,k,URM_Tr_,U_old_,H_old_, V_old_) ;				
			for (uword y = 0 ; y < U_.n_cols ; ++y)
			{
	#ifdef AMFOPENMP
	#pragma omp atomic
	#endif
				G(x,y) = G(x,y) - aux * A(y,k);
			}
		}

	#ifdef AMFTIME
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << "ms" << std::endl;
	#endif

	// Run the loop

	for (unsigned int n=0 ; n < n_max_iter_gradient_  ; ++n ) 
	{
		solve_pg_U_One_Iteration(G,A,AAt);
	}

	#ifndef NDEBUG
	double obj = evaluate_Obj_Function(URM_Tr_,U_,H_old_,V_old_,U_old_,H_old_,V_old_,lambda_);
	std::cout << "Objective function : " << obj << std::endl;
	#endif

}

void AMF::solve_pg_U_With_Log()
{
	std::ofstream logfile;
	logfile.open(amf_filename_prefix_+"log_pg_u.txt");

	// Initialize U (with "warm-up")
	U_ = U_old_;

	mat A = H_old_*V_old_;
	mat AAt = A*(A.t() );

	mat G(U_.n_rows,U_.n_cols,fill::zeros);

	gradient_step_ = 1.0;

	// Compute the constant part of the gradient
	#ifndef NDEBUG
	std::cout << "Computing linear part of gradient ... " ;
	#endif

	#ifdef AMFTIME 
	auto begin = std::chrono::high_resolution_clock::now();
	#endif

	#ifdef AMFOPENMP
	#pragma omp parallel for
	#endif
	for (uword k = 0 ; k < A.n_cols ; ++k)
		for (uword x =0 ; x< U_.n_rows ; ++x)
		{
			double aux = build_S(x,k,URM_Tr_,U_old_,H_old_, V_old_) ;				
			for (uword y = 0 ; y < U_.n_cols ; ++y)
			{
	#ifdef AMFOPENMP
	#pragma omp atomic
	#endif
				G(x,y) = G(x,y) - aux * A(y,k);
			}
		}

	#ifdef AMFTIME
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << "ms" << std::endl;
	#endif

	// Run the loop

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
	// Update gradient with non-constant part
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

	#ifndef NDEBUG
	std::cout << "Step = " << gradient_step_ << " Value = " << res << " ";
	std::cout <<( (res <= 0)?("First is feasible so I increase the step"):("First is not feasible") )<< std::endl;
	#endif

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

		#ifndef NDEBUG
		std::cout << "Step = " << gradient_step_ << " Value = " << res << " ";
		std::cout <<( (res <= 0)?("Feasible"):("Not feasible") )<< std::endl;
		#endif


		if (res <=0)
			is_feasible = true;
		else 
			gradient_step_ = beta*gradient_step_;
	}

	#ifndef NDEBUG
	std::cout << "Selected step is " << gradient_step_ << std::endl;
	#endif

	// Update U_ 
	U_ = U_ - gradient_step_*G -gradient_step_*lambda_*U_ ;
	get_Positive_Matrix(U_);


}

void AMF::solve_pg_H()
{
    H_ = H_old_;

    mat UtU = (U_.t())*U_;

	#ifdef AMFTIME 
	std::cout << "Computing V*Vt ... " ;
	auto begin_VVt = std::chrono::high_resolution_clock::now();
	#endif

    mat VVt(V_old_.n_rows,V_old_.n_rows,fill::zeros);
	for (uword i = 0 ; i < V_old_.n_rows ; ++i )
		for (uword j = 0 ; j < V_old_.n_rows ; ++j )
		{
			for (uword k = 0 ; k < V_old_.n_cols ; ++k )
				VVt(i,j) += V_old_(i,k)*V_old_(j,k);
		}


	#ifdef AMFTIME
	auto end_VVt = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_VVt-begin_VVt).count() << "ms" << std::endl;
	#endif

	gradient_step_ = 1.0;

    mat G(H_.n_rows,H_.n_cols,fill::zeros);

    #ifndef NDEBUG
    std::cout << "Computing linear part of gradient ..." << std::endl;
    #endif

	#ifdef AMFTIME 
	auto begin_LH = std::chrono::high_resolution_clock::now();
	#endif


 	#ifdef AMFOPENMP
 	#pragma omp parallel for
 	#endif
	for (uword k = 0 ; k < U_.n_rows ; ++k)
 	for (uword l = 0 ; l < V_old_.n_cols ; ++l)
 	{
 		double aux = build_S(k,l,URM_Tr_,U_old_,H_old_, V_old_);
		for (uword x =0 ; x< H_.n_rows ; ++x)
		for (uword y = 0 ; y < H_.n_cols ; ++y)
		{
	#ifdef AMFOPENMP
 	#pragma omp atomic
 	#endif
 				G(x,y) = G(x,y) - U_(k,x)*aux*V_old_(y,l) ;
 		}
       	
    }

    #ifdef AMFTIME
	auto end_LH = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_LH-begin_LH).count() << "ms" << std::endl;
	#endif
    
    for (unsigned int n=0 ; n < n_max_iter_gradient_ ; ++n )
    {
        solve_pg_H_One_Iteration(G,UtU,VVt);
    }

    #ifndef NDEBUG
    double obj = evaluate_Obj_Function(URM_Tr_,U_,H_,V_old_,U_old_,H_old_,V_old_,lambda_);
    std::cout << "Objective function : " << obj << std::endl;
    #endif

}

void AMF::solve_pg_H_With_Log(){

    std::ofstream logfile(amf_filename_prefix_+"log_pg_h.txt");

    // Initialize H (with warm up)
    H_ = H_old_;

    // Precompute UtU and VVt

    mat UtU = (U_.t())*U_;


	#ifdef AMFTIME 
	std::cout << "Computing V*Vt ... " ;
	auto begin_VVt = std::chrono::high_resolution_clock::now();
	#endif

    mat VVt(V_old_.n_rows,V_old_.n_rows,fill::zeros);
	for (uword i = 0 ; i < V_old_.n_rows ; ++i )
		for (uword j = 0 ; j < V_old_.n_rows ; ++j )
		{
			for (uword k = 0 ; k < V_old_.n_cols ; ++k )
				VVt(i,j) += V_old_(i,k)*V_old_(j,k);
		}


	#ifdef AMFTIME
	auto end_VVt = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_VVt-begin_VVt).count() << "ms" << std::endl;
	#endif

	// Initialize gradient step
	gradient_step_ = 1.0;

    // Precompute constant part of the gradient
    mat G(H_.n_rows,H_.n_cols,fill::zeros);

    #ifndef NDEBUG
    std::cout << "Computing linear part of gradient ..." << std::endl;
    #endif

	#ifdef AMFTIME 
	auto begin_LH = std::chrono::high_resolution_clock::now();
	#endif

 	#ifdef AMFOPENMP
 	#pragma omp parallel for
 	#endif
	for (uword k = 0 ; k < U_.n_rows ; ++k)
 	for (uword l = 0 ; l < V_old_.n_cols ; ++l)
 	{
 		double aux = build_S(k,l,URM_Tr_,U_old_,H_old_, V_old_);

		for (uword x =0 ; x< H_.n_rows ; ++x)
		for (uword y = 0 ; y < H_.n_cols ; ++y)
		{
	#ifdef AMFOPENMP
	#pragma omp atomic
	#endif
 				G(x,y) = G(x,y) - U_(k,x)*aux*V_old_(y,l) ;
 		}
       	
    }

    #ifdef AMFTIME
	auto end_LH = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_LH-begin_LH).count() << "ms" << std::endl;
	#endif
    
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

	#ifndef NDEBUG
	std::cout << "Step = " << gradient_step_ << " Value = " << res << " ";
	std::cout <<( (res <= 0)?("First is feasible so I increase the step"):("First is not feasible") )<< std::endl;
	#endif

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

		#ifndef NDEBUG
		std::cout << "Step = " << gradient_step_ << " Value = " << res << " ";
		std::cout <<( (res <= 0)?("Feasible"):("Not feasible") )<< std::endl;
		#endif

		if (res <=0)
			is_feasible = true;
		else 
			gradient_step_ = beta*gradient_step_;
	}

	#ifndef NDEBUG
	std::cout << "Selected step is " << gradient_step_ << std::endl;
	#endif

    // Update H_
    H_ = H_ - gradient_step_*G;
    get_Positive_Matrix(H_);

}

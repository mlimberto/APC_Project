//
//  compute_V.cpp
//
//
//
//  Authors : Lara Bombardieri, Andrea Brandoli, Matteo Limberto

#include "amf.h"

using namespace arma;

//////////////////////////
////  METHODS       //////
//////////////////////////


void AMF::solve_V()
{

    V_ = V_old_;

    #ifndef NDEBUG
    std::cout << "Computing W = UH ..." << std::endl;
    #endif

    mat W = U_*H_;
    
    #ifndef NDEBUG
    std::cout << "Computing WtW" << std::endl;
    #endif

    mat WtW = (W.t())*W;

    mat G(V_.n_rows,V_.n_cols,fill::zeros);

    gradient_step_ = 1e-4;

    // Compute the constant part of the gradient
    #ifndef NDEBUG
    std::cout << "Computing linear part of gradient..." << std::endl;
    #endif

    #ifdef AMFTIME 
    auto begin_LV = std::chrono::high_resolution_clock::now();
    #endif

    #ifdef AMFOPENMP
    #pragma omp parallel for 
    #endif
    for (uword y = 0; y < G.n_cols ; ++y)
    {
        for (uword k=0 ; k< W.n_rows ; ++k)
        {
            double aux = build_S(k,y,URM_Tr_,U_old_,H_old_, V_old_);
            for (uword x = 0; x < G.n_rows ; ++x)
            {
                G(x,y) = G(x,y) - W(k,x)*aux;
            }
        }
    }

    #ifdef AMFTIME
    auto end_LV = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_LV-begin_LV).count() << "ms" << std::endl;
    #endif

    // Run the loop

    bool stop_criterion = false;

    double prec_obj(0) , curr_obj(0);

    curr_obj = evaluate_Obj_Function(URM_Tr_,U_,H_,V_old_,U_old_,H_old_,V_old_,lambda_);

    #ifndef NDEBUG
    std::cout << "Evaluating first objective function : " << curr_obj << std::endl;
    #endif

    for (unsigned int n=0 ; (n < n_max_iter_gradient_ ) && (!stop_criterion) ; ++n ) 
    {
        prec_obj = curr_obj;
        solve_V_One_Iteration(G,WtW,prec_obj,curr_obj);

        // Evaluate stop criterion

        if (abs(curr_obj - prec_obj)/curr_obj < toll_gradient_)
            stop_criterion = true;

        #ifndef NDEBUG
        std::cout << "old_obj " << prec_obj << " new_obj" << curr_obj << std::endl;
        #endif

        // Print information
        std::cout << "Iteration " << n+1 << " : Objective function = " << curr_obj << std::endl;

    }

}

void AMF::solve_V_With_Log()
{
    std::ofstream logfile;
    logfile.open(amf_filename_prefix_+"log_pg_v.txt");

    // Initialize V (with "warm-up")
    V_ = V_old_;

    // Compute W and WtW

    #ifndef NDEBUG
    std::cout << "Computing W = UH ..." << std::endl;
    #endif

    mat W = U_*H_;
    
    #ifndef NDEBUG
    std::cout << "Computing WtW" << std::endl;
    #endif

    mat WtW = (W.t())*W;

    // Initialize the gradient matrix and set the gradient step

    mat G(V_.n_rows,V_.n_cols,fill::zeros);

    gradient_step_ = 1e-4;

    // Compute the constant part of the gradient
    #ifndef NDEBUG
    std::cout << "Computing linear part of gradient..." << std::endl;
    #endif

    #ifdef AMFTIME 
    auto begin_LV = std::chrono::high_resolution_clock::now();
    #endif

    #ifdef AMFOPENMP
    #pragma omp parallel for 
    #endif
    for (uword y = 0; y < G.n_cols ; ++y)
    {
        for (uword k=0 ; k< W.n_rows ; ++k)
        {
            double aux = build_S(k,y,URM_Tr_,U_old_,H_old_, V_old_);
            for (uword x = 0; x < G.n_rows ; ++x)
            {
                G(x,y) = G(x,y) - W(k,x)*aux;
            }
        }
    }

    #ifdef AMFTIME
    auto end_LV = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_LV-begin_LV).count() << "ms" << std::endl;
    #endif

    // Run the loop

    bool stop_criterion = false;

    double prec_obj(0) , curr_obj(0);

    curr_obj = evaluate_Obj_Function(URM_Tr_,U_,H_,V_old_,U_old_,H_old_,V_old_,lambda_);

    #ifndef NDEBUG
    std::cout << "Evaluating first objective function : " << curr_obj << std::endl;
    #endif

    for (unsigned int n=0 ; (n < n_max_iter_gradient_ ) && (!stop_criterion) ; ++n ) 
    {
        prec_obj = curr_obj;
        solve_V_One_Iteration(G,WtW,prec_obj,curr_obj);

        // Evaluate stop criterion

        if (abs(curr_obj - prec_obj)/curr_obj < toll_gradient_)
            stop_criterion = true;

        #ifndef NDEBUG
        std::cout << "old_obj " << prec_obj << " new_obj " << curr_obj << std::endl;
        #endif

        // Save information on logfile
        logfile << curr_obj << "\n";
        total_logfile_ << curr_obj << "\n";

        // Print information
        std::cout << "Iteration " << n+1 << " : Objective function = " << curr_obj << std::endl;

    }


    logfile.close();

}

void AMF::solve_V_One_Iteration(arma::mat G,const arma::mat &WtW, const double prec_obj,double &curr_obj)
{
    // Update gradient with non-constant part
    G += WtW*V_;

    bool is_feasible = false;
    mat V_hat = V_ - gradient_step_*G ;
    sp_mat V_temp=V_;
    project_V(V_temp,V_hat);

    curr_obj = evaluate_Obj_Function(URM_Tr_,U_,H_,V_temp,U_old_,H_old_,V_old_,lambda_);

    //Check if the object function increases
    if (curr_obj < prec_obj){
        // Update V
        // V_=V_temp;
        // is_feasible=true;
        #ifndef NDEBUG
        std::cout << gradient_step_ << " is feasible, increasing gradient step" << std::endl;
        #endif
        gradient_step_ = gradient_step_/0.1;
    }
    else{
        #ifndef NDEBUG
        std::cout << gradient_step_ << " not feasible, Decreasing gradient step" << std::endl;
        #endif
        gradient_step_ = gradient_step_ *0.1;
    }


    while(!is_feasible){
        V_hat = V_ - gradient_step_*G ;
        V_temp=V_;
        project_V(V_temp,V_hat);

        curr_obj = evaluate_Obj_Function(URM_Tr_,U_,H_,V_temp,U_old_,H_old_,V_old_,lambda_);

        if (curr_obj < prec_obj){
            // Update V
            V_=V_temp;
            is_feasible=true;
            #ifndef NDEBUG
            std::cout << gradient_step_ << " is feasible" << std::endl;
            #endif

        }else{
            #ifndef NDEBUG
            std::cout <<gradient_step_ << "not feasible, Decreasing gradient step" << std::endl;
            #endif
            gradient_step_ = gradient_step_ *0.1;
        }

        if (gradient_step_ < 1e-10)
            break;

    }

}     

void AMF::project_V(arma::sp_mat &V_new ,arma::mat &V_hat)
{
    // Projection on ICM
    sp_mat T_hat=project_ICM(V_hat);


    // Projection on SUM = 1
    uvec sorted_indices;

    for(uword j=0; j<T_hat.n_cols;++j){
        uword n_nonzero=T_hat.col(j).n_nonzero;
        sorted_indices=sort_index(mat(T_hat).col(j),"descend");
        vec t_hat=sort(mat(T_hat).col(j),"descend");

        double tau=0;
        double tau_old=0;
        for(uword i=0;i<n_nonzero;++i){
            // Computing tau
            tau_old=(double)(sum(t_hat(span(0,i)))-1)/(i+1);
            if(tau_old<t_hat(i)){
                tau=tau_old;
            }
        }
        uword i=0;
        while(i<n_nonzero && t_hat(i)-tau>0){
             V_new(sorted_indices(i),j)=t_hat(i)-tau;
            i++;
        }
    }

}



sp_mat AMF::project_ICM(const mat &G){

    // Funtion to set to zero the elements of the matrix B where ICM is zero

    if (G.n_cols != ICM_.n_cols || G.n_rows != ICM_.n_rows)    {
        std::cerr << "ERROR in project_ICM : Inconsistent matrices !!!" << std::endl;
    }

    sp_mat P(ICM_Location_Matrix_,G.elem(get_Vector_Of_Indices(ICM_Location_Matrix_)));
    return P;
}



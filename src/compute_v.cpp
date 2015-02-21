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

    // Compute the linear part of the gradient
    #ifndef NDEBUG
    std::cout << "Computing linear part of gradient..." << std::endl;
    #endif

    #ifdef AMFTIME 
    auto begin_LV = std::chrono::high_resolution_clock::now();
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

    // Compute the linear part of the gradient
    #ifndef NDEBUG
    std::cout << "Computing linear part of gradient..." << std::endl;
    #endif

    #ifdef AMFTIME 
    auto begin_LV = std::chrono::high_resolution_clock::now();
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
    // Update gradient with quadratic part
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
        std::cout << gradient_step_ << " is feasible, increasing gradient step" << std::endl;
        gradient_step_ = gradient_step_/0.1;
    }
    else{
        std::cout << gradient_step_ << " not feasible, Decreasing gradient step" << std::endl;
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
            std::cout << gradient_step_ << " is feasible" << std::endl;

        }else{
            std::cout <<gradient_step_ << "not feasible, Decreasing gradient step" << std::endl;
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




sp_mat AMF::solve_V_One_Step_Gradient(const sp_mat &V_0){

    std::cout<<"Solving One gradient Step"<<std::endl;

    // Computation of the gradient of f=||S-UHV||^2
    std::cout<<"Solving One gradient Step : computing gradient"<<std::endl;
    mat UH = U_*H_;
    mat grad=2*UH.t()*(UH*V_0-project_URM(URM_Tr_,U_old_*H_old_*V_old_));

    // Projection of the gradient on ICM
    std::cout<<"Solving One gradient Step : projecting on ICM"<<std::endl;
    sp_mat grad_tilda=project_ICM(grad);
   // grad_tilda.print("grad = ");

    // Projection of the gradient on the matrix space whose columns sum to zero
    std::cout<<"Solving One gradient Step : orthogonal_projection"<<std::endl;
    orthogonal_projection(grad_tilda);
    //grad_tilda.print("grad = ");

    // Computation of t
    std::cout<<"Solving One gradient Step : computing t"<<std::endl;
    double t=1/norm(UH.t()*UH,2);
    for(uword alpha=0 ; alpha<grad_tilda.n_rows ; alpha++){
        for(uword beta=0 ; beta<grad_tilda.n_cols ; beta++){
            if(grad_tilda(alpha,beta)>0){
                if (V_0(alpha,beta)==0){
                    std::cerr<<"ERROR : V_0("<<alpha<<","<<beta<<")=0 but grad_tilda("<<alpha<<","<<beta<<")>0 !!"<<std::endl;
                }else{
                    t=std::min(t,V_0(alpha,beta)/grad_tilda(alpha,beta));
                }
            }
        }
    }
    std::cout<<"t = "<<t<<std::endl;
    // Computation of V^
    return V_0-gradient_step_*grad_tilda;

}

sp_mat AMF::project_ICM(const mat &G){

    // Funtion to set to zero the elements of the matrix B where ICM is zero

    if (G.n_cols != ICM_.n_cols || G.n_rows != ICM_.n_rows)    {
        std::cerr << "ERROR in project_ICM : Inconsistent matrices !!!" << std::endl;
    }

    sp_mat P(ICM_Location_Matrix_,G.elem(get_Vector_Of_Indices(ICM_Location_Matrix_)));
    return P;
}

void AMF::orthogonal_projection(sp_mat &G){

    for(uword j=0 ; j<G.n_cols ; j++)
    {

        auto v = G.col(j);    

        sp_mat u(v.n_rows,v.n_cols);
        uword n_nonzero=G.col(j).n_nonzero;
        for (uword i = 0 ; i<v.n_rows ; i++){
            if (v(i,0)!=0){
       //         n_nonzero+=1;
              u(i,0)=1;
            }
        }

        if (n_nonzero!=0)
        {
            u=u/(sqrt(n_nonzero));
            v=v-as_scalar(v.t()*u)*u;
        }
    }

}

sp_mat AMF::solve_V_One_Step_Gradient2(const sp_mat &V_0){
    std::cout<<"Solving One gradient Step2"<<std::endl;

    // Computation of the gradient of f=||S-UHV||^2
    std::cout<<"Solving One gradient Step2 : computing gradient"<<std::endl;
    mat UH = U_*H_;
    mat grad=2*UH.t()*(UH*V_0-project_URM(URM_Tr_,U_old_*H_old_*V_old_));

    //When dimensions are too big, do it by column :
    //mat grad(V_old_.n_rows,V_old_.n_cols);
    //for(uword j=0; j<V_old_.n_cols ; ++j){
    //    grad.col(j)=2*UH.t()*(UH*V_old_.col(j)-project_URM_by_column(j,U_old_*H_old_*V_old_.col(j)));
    //}

    double t=1/norm(UH.t()*UH,2);
    mat V_hat=V_0-t*grad;

    // Projection on ICM
    std::cout<<"Solving One gradient Step : projecting V_hat on ICM"<<std::endl;
    sp_mat T_hat=project_ICM(V_hat);
    sp_mat V_new(V_0.n_rows,V_0.n_cols);

    uvec sorted_indices;
    std::cout<<"Solving One gradient Step : projecting on S"<<std::endl;

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
    return V_new;

}

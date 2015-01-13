//
//  compute_V.cpp
//
//
//
//

#include "amf.h"


//////////////////////////
////  METHODS       //////
//////////////////////////


void AMF::solve_V(){

    // Computation of V
    #ifndef NDEBUG
    std::cout << "Computing V " << std::endl;
    #endif
    sp_mat V_hat=solve_V_One_Step_Gradient(V_old_);

}

sp_mat AMF::solve_V_One_Step_Gradient(const sp_mat &V_0){

    std::cout<<"Solving One gradient Step"<<std::endl;

    // Computation of the gradient of f=||S-UHV||^2
    std::cout<<"Solving One gradient Step : computing gradient"<<std::endl;
    mat UH = U_*H_;
    mat grad=2*UH.t()*(UH*V_0-project_URM(URM_,U_old_*H_old_*V_old_));

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
    return V_0-t*grad_tilda;

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
    mat grad=2*UH.t()*(UH*V_0-project_URM(URM_,U_old_*H_old_*V_old_));

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

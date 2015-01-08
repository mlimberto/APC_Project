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

    // Computation of the gradient of f=||S-UHV||^2

    mat UH = U_*H_;
    mat grad(V_0.n_rows,V_0.n_cols,fill::zeros);
    double sum(0);
    for (uword alpha = 0 ; alpha<grad.n_rows ; alpha++){
        for(uword beta=0 ; beta<grad.n_cols ; beta++){
            sum=0;
            for(uword i=0 ; i<UH.n_rows ; i++){
                sum+=UH(i,alpha)*(UH.row(i)*V_0.col(beta)-build_S(i,beta,URM_,U_old,H_old,V_old));
            }
            grad(alpha,beta)=2*sum;
        }
    }

    // Projection of the gradient on C
    sp_mat grad_tilda=project_ICM(grad);

    // Projection of the gradient on the matrix space whose columns sum to zero
    for(uword j=0 ; j<grad_tilda.cols ; j++){
        orthogonal_projection(grad_tilda.col(j));
    }

    // Computation of t
    double t=1/norm(UH.t()*UH,2);
    for(uword alpha=0 ; alpha<grad_tilda.rows ; alpha++){
        for(uword beta=0 ; beta<grad_tilda.cols ; beta++){
            if(grad_tilda(alpha,beta)>0){
                if (V_0(alpha,beta)==0){
                    std::cerr<<"ERROR : V_0("<<alpha<<","<<beta<<")=0 but grad_tilda("<<alpha<<","<<beta<<")>0 !!"<<std::endl;
                }else{
                    t=min(t,V_0(alpha,beta)/grad_tilda(alpha,beta));
                }
            }
        }
    }
    // Computation of V^
    return V_0-t*grad_tilda;

}

sp_mat AMF::project_ICM(const mat &G){

    // Funtion to set to zero the elements of the matrix B where ICM is zero

    if (G.n_cols != ICM_.n_cols || G.n_rows != ICM_.n_rows)
    {
        std::cerr << "ERROR in project_ICM : Inconsistent matrices !!!" << std::endl;
    }

    sp_mat A(ICM_.n_rows,ICM_.n_cols);

    for (uword i = 0 ; i<ICM_.n_rows ; ++i){
        for (uword j = 0 ; j < ICM_.n_cols ; ++j){
            A(i,j)=G(i,j)*(ICM_(i,j)>0);
        }
    }
    return A;
}

void AMF::othogonal_projection(sp_mat &v){

    // Check that v is a column vector
    if(!v.is_colvec()){
        std::cerr << "ERROR in othogonal_projection : v must be a column vector !!!" << std::endl;
    }
    sp_mat u(v.rows,v.cols,fill::zeros);
    uword n_nonzero(0);
    for (uword i = 0 ; i<v.n_rows ; i++){
        if (v(i,0)!=0){
            n_nonzero+=1;
            u(i,0)=1;
        }
    }
    u=u/(sqrt(n_nonzero));
    v=v-as_scalar(v*u)*u;
}

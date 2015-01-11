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
    mat grad(V_0.n_rows,V_0.n_cols,fill::zeros);
    //double sum(0);
    //for (uword alpha = 0 ; alpha<grad.n_rows ; alpha++){
      //  for(uword beta=0 ; beta<grad.n_cols ; beta++){
        //    sum=0;
          //  for(uword i=0 ; i<UH.n_rows ; i++){
            //    sum+=UH(i,alpha)*(as_scalar(UH.row(i)*V_0.col(beta))-build_S(i,beta,URM_,U_old_,H_old_,V_old_));
           // }
            //grad(alpha,beta)=2*sum;
       // }
    //}

    //for (uword alpha = 0 ; alpha<grad.n_rows ; alpha++){
      //  for(uword beta=0 ; beta<grad.n_cols ; beta++){
        //    std::cout<<"("<<alpha<<","<<beta<<")"<<std::endl;
          //  grad(alpha,beta)=2*dot(UH.col(alpha),(UH*V_0.col(beta)-build_S_by_column(beta,URM_,U_old_,H_old_,V_old_)));
        //}
    //}

    grad=2*UH.t()*(UH*V_0-project_URM(URM_,U_old_*H_old_*V_old_));
    //for (uword alpha = 0 ; alpha<grad.n_rows ; alpha++){
      //  std::cout<<alpha<<std::endl;
        //grad.row(alpha)=UH.col(alpha).t()*(UH*V_0-S);
    //}

    // Projection of the gradient on C
    std::cout<<"Solving One gradient Step : projecting on C"<<std::endl;
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

void AMF::orthogonal_projection(sp_mat &G){

    for(uword j=0 ; j<G.n_cols ; j++)
    {

        auto v = G.col(j);    

        sp_mat u(v.n_rows,v.n_cols);
        uword n_nonzero(0);
        for (uword i = 0 ; i<v.n_rows ; i++){
            if (v(i,0)!=0){
                n_nonzero+=1;
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
    mat grad(V_0.n_rows,V_0.n_cols,fill::zeros);
    mat V_hat=V_0-gradient_step_*grad;
    sp_mat T_hat=project_ICM(V_hat);
    sp_mat t_hat;
    uvec sorted_indeces;
    for(uword j=0; j<T_hat.n_cols;++j){
        sorted_indeces=sort_index(T_hat.col(j),"descend");
        t_hat=sort(T_hat.col(j),"descend");
        uword i=0;
        double tau=0;
        doubel tau_old=0;
        while(t_hat(i)!=0){
            tau_old=t_hat()
            i++;
        }
    }

}

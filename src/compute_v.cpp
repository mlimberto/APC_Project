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

}

void AMF::solve_V_One_Step_Gradient(){

    // Computation of the gradient of f=||S-UHV||^2

    mat UH = U_*H_;

    sp_mat grad(V_old_.n_rows,V_old_.n_cols); // Abbiamo automaticamente l'inizializzazione a zero


    // Projection of the gradient on C

    // Projection of the gradient on the matrix space whose columns sum to zero

    // Computation of t

    // Computation of V^
}

sp_mat AMF::project_ICM(const mat &B){

    if (B.n_cols != ICM_.n_cols || B.n_rows != ICM_.n_rows)
    {
        std::cerr << "ERROR in project_ICM : Inconsistent matrices !!!" << std::endl;
    }
    sp_mat A(ICM_.n_cols,ICM_.n_rows);

    for (auto i = 0 ; i<ICM_.n_rows ; ++i){
        for (auto j = 0 ; j < ICM_.n_cols ; ++j){
            A(i,j)=B(i,j)*(ICM_(i,j)>0);
        }
    }
    return A;
}

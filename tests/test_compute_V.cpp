#include <iostream>
#include "amf.h"

#include <armadillo>

using namespace std;
using namespace arma;

bool test_orthogonal_projection()
{
    //sp_mat A= sprandu<sp_mat>(10,15,0.1);
    // A(0,0) = -2;
    //A(1,0) = 1;
    //A.print("Matrix A : ");

    AMF amf;

    //initialize URM and ICM randomly (just for test)

    //amf.initialize_Parameters("test_parameters.txt") ;
    amf.set_lambda(1.5);
    amf.set_n_max_iter(10000);
    amf.set_toll(1e-3);
    amf.set_n_latent_factors(50);
    amf.set_n_max_iter_gradient(100);
    amf.set_toll_gradient(1e-2);
    amf.set_gradient_step(0.1);

    cout << "Lambda has been set to " << amf.get_lambda() << endl;
    cout << "Lat factors has been set to " << amf.get_n_latent_factors() << endl;
    cout << "Toll has been set to " << amf.get_toll() << endl;
    cout << "n_max_iter has been set to " << amf.get_n_max_iter() << endl;
    cout << "n_max_iter_gradient has been set to " << amf.get_n_max_iter_gradient() << endl;
    cout << "toll_gradient has been set to " << amf.get_toll_gradient() << endl;
    cout << "gradient_step has been set to " << amf.get_gradient_step() << endl;
    cout << ( (amf.get_lambda()==1.5) ? "Succeed" : "Failed") << endl;

    amf.initialize_URM_Locations("../APC_Project/dataset/Sampled_Dataset/urm_sampling.txt");
    //amf.print_URM();
    amf.initialize_ICM_Locations("../APC_Project/dataset/Sampled_Dataset/icm_sampling.txt");
    //amf.print_ICM();
    cout<<"URM size = "<<amf.URM_.n_rows<<" x "<<amf.URM_.n_cols<<endl;
    cout<<"ICM size = "<<amf.ICM_.n_rows<<" x "<<amf.ICM_.n_cols<<endl;

    //Inizialize matrices (randomly)
    amf.U_=randu<mat>(amf.URM_.n_rows,amf.get_n_latent_factors());
    amf.H_=randu<mat>(amf.get_n_latent_factors(),amf.ICM_.n_rows);
    amf.U_old_=randu<mat>(amf.URM_.n_rows,amf.get_n_latent_factors());
    amf.H_old_=randu<mat>(amf.get_n_latent_factors(),amf.ICM_.n_rows);
    amf.V_old_=speye<sp_mat>(amf.ICM_.n_rows,amf.ICM_.n_cols);
    cout<<"object function before gradient = "<<evaluate_Obj_Function(amf.URM_,amf.U_,amf.H_,amf.V_old_,amf.U_old_,amf.H_old_,amf.V_old_,amf.lambda_)<<endl;
    sp_mat V_hat=amf.solve_V_One_Step_Gradient(amf.V_old_);
    cout<<"object function after gradient = "<<evaluate_Obj_Function(amf.URM_,amf.U_,amf.H_,V_hat,amf.U_old_,amf.H_old_,amf.V_old_,amf.lambda_)<<endl;
    sp_mat diff=amf.V_old_-V_hat;
    diff.print("V_old_-V_hat : ");
    //amf.orthogonal_projection(A);
    //A.print("Resulting matrix");



    return true;
}

int main(int argc,char **argv)
{
    mat B = randu<mat>(4,5);
    B.print("Matrix B : ");
    sp_mat ICM_example=sprandu<sp_mat>(B.n_rows,B.n_cols, 0.1);
    ICM_example.print("ICM_example : ");

    sp_mat A(B.n_rows,B.n_cols);

    for (uword i = 0 ; i<ICM_example.n_rows ; ++i){
        for (uword j = 0 ; j < ICM_example.n_cols ; ++j){
            A(i,j)=B(i,j)*(ICM_example(i,j)>0);
        }
    }
    A.print("A (projection of B on ICM_example) : ");
    //sp_mat D=A/2;
    //D.print("D :");
    //sp_mat A2(B.n_rows,B.n_cols);
    //A2=B(i,j)%(ICM_example>0);
    //A2.print("A2 (projection of B on ICM_example) : ");
    //sp_mat v=A.col(0);
    //v.print("v : ");

    // Test per la proiezione ortogonale
    test_orthogonal_projection();


    return 0;
}

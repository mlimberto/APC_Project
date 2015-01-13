#include <iostream>
#include "amf.h"

#include <armadillo>

using namespace std;
using namespace arma;

bool test_orthogonal_projection(){


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
    //amf.U_old_ = mat(amf.URM_.n_rows,amf.get_n_latent_factors(),fill::ones);
    amf.H_old_=randu<mat>(amf.get_n_latent_factors(),amf.ICM_.n_rows);
    //amf.H_old_ = mat(amf.get_n_latent_factors(),amf.ICM_.n_rows,fill::eye);
    amf.initialize_matrices();
    cout<<"Checking constraiV_old_nts on V_old_..."<<endl;
    check_V_Constraint(amf.V_old_,amf.ICM_);


    cout<<"Evaluating object function ..."<<endl;
    double before=evaluate_Obj_Function(amf.URM_,amf.U_,amf.H_,amf.V_old_,amf.U_old_,amf.H_old_,amf.V_old_,amf.lambda_);
    cout<<" Object function before gradient = "<<before<<endl;

    sp_mat V_new1=amf.solve_V_One_Step_Gradient(amf.V_old_);

    cout<<"Evaluating object function ..."<<endl;
    cout<<"Object function after gradient = "<<evaluate_Obj_Function(amf.URM_,amf.U_,amf.H_,V_new1,amf.U_old_,amf.H_old_,amf.V_old_,amf.lambda_)<<endl;

    cout<<"Checking constraints on V_new1..."<<endl;
    check_V_Constraint(V_new1,amf.ICM_);

    sp_mat V_new2=amf.solve_V_One_Step_Gradient2(amf.V_old_);

    cout<<"Evaluating object function ..."<<endl;
    cout<<"Object function after gradient = "<<evaluate_Obj_Function(amf.URM_,amf.U_,amf.H_,V_new2,amf.U_old_,amf.H_old_,amf.V_old_,amf.lambda_)<<endl;

    cout<<"Checking constraints on V_new2..."<<endl;
    check_V_Constraint(V_new2,amf.ICM_);




    return true;
}


int main(int argc,char **argv)
{

    umat locations;
    locations << 0 << 1 << 3 << endr
              << 0 << 2 << 4 << endr;

    vec values;
    values << 1 << 5 << 1 << endr;

    sp_mat ICM_example(locations, values);
    ICM_example.print("ICM_example : ");
    mat B = randu<mat>(ICM_example.n_rows,ICM_example.n_cols);
    B.print("Matrix B : ");
    uvec q1 = find(mat(ICM_example));
    q1.print("q1 :");
    uvec v=get_Vector_Of_Indices(locations);
    v.print("v : ");
    sp_mat P(locations,B.elem(get_Vector_Of_Indices(locations)));
    P.print("B projected on ICM: = ");
    //cout<<"number of non zero = "<<P.n_nonzero<<endl;
    // Test per la proiezione ortogonale
    test_orthogonal_projection();


    return 0;
}


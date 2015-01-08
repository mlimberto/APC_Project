#include <iostream>

#include <armadillo>

using namespace std;
using namespace arma;

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
    sp_mat D=A/2;
    D.print("D :");
    //sp_mat A2(B.n_rows,B.n_cols);
    //A2=B(i,j)%(ICM_example>0);
    //A2.print("A2 (projection of B on ICM_example) : ");
    sp_mat v=A.col(0);
    v.print("v : ");
    return 0;
}

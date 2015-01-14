#include <iostream> 
#include <string>
#include <vector>

#include "my_utils.h"

using namespace std;
using namespace arma;

int main(int argc,char** argv)
{

    std::string matrixICM_filename("../dataset/Sampled_Dataset/icm_sampling.txt");

    std::string matrixURM_filename("../dataset/Sampled_Dataset/urm_sampling.txt");

    sp_umat M;
    umat location_matrix;
    Col<uword> values;

    import_Sparse_Matrix<uword>(matrixICM_filename,M,location_matrix,values);

    M.print();

    SpMat<double> urm;
    Col<double> val;

    import_Sparse_Matrix<double>(matrixURM_filename,urm,location_matrix,val);

    urm.print();


    return 0;
}

//
//  my_utils.cpp
//  
//
//
//

#include "my_utils.h"

using std::string;
using std::vector;

using namespace arma;


// ricordarsi che gli indici di riga e colonna in c++ partono da 0 !!

template<typename T>
bool import_Sparse_Matrix(std::string mfilename,SpMat<T> &MM,umat &location_mat,Col<T> &values)
{
    arma::umat RCi(2,1);
    T val;
    std::ifstream matrix_file(mfilename);

    location_mat.reset();
    values.reset();
    MM.reset();

    // Read the file and build the Location Matrix and the Values vector
    unsigned int i=0;
    while (matrix_file >> RCi(0,0) >> RCi(1,0) >> val){
        location_mat.insert_cols(i,RCi);
        values.resize(i+1);
        values(i)=val;
        i++;
    }

    MM = SpMat<T>(location_mat,values);

	return true;
}
template bool import_Sparse_Matrix<uword>(std::string mfilename,SpMat<uword> &MM,umat &location_mat,Col<uword> &values);
template bool import_Sparse_Matrix<double>(std::string mfilename,SpMat<double> &MM,umat &location_mat,Col<double> &values);


template<typename T>
bool import_Sparse_Matrix(std::string mfilename,SpMat<T> &MM)
{
    arma::umat RCi(2,1);
    T val;
    std::ifstream matrix_file(mfilename);

    arma::umat location_mat;
    Col<T> values;
    MM.reset();

    // Read the file and build the Location Matrix and the Values vector
    unsigned int i=0;
    while (matrix_file >> RCi(0,0) >> RCi(1,0) >> val){
        location_mat.insert_cols(i,RCi);
        values.resize(i+1);
        values(i)=val;
        i++;
    }

    MM = SpMat<T>(location_mat,values);

	return true;
}
template bool import_Sparse_Matrix<uword>(std::string mfilename,SpMat<uword> &MM);
template bool import_Sparse_Matrix<double>(std::string mfilename,SpMat<double> &MM);



double build_S(uword i, uword j,const SpMat<double>& URM,const Mat<double>& U,const Mat<double>& H,const SpMat<double>& V )
{
	if (URM(i,j) != 0){
		return URM(i,j);
	}
	else
	{
		return as_scalar(U.row(i)*H*V.col(j));
	}
}


mat build_S_by_column(uword j,const SpMat<double>& URM,const Mat<double>& U,const Mat<double>& H,const SpMat<double>& V ){
    mat c=(U*H*V.col(j));
    for(uword i=0; i<URM.n_rows; ++i){
        if (URM(i,j) != 0){
            c(i,0)=URM(i,j);
        }
    }
    return c;
}


mat project_URM(const SpMat<double>& URM, const mat &S){
    mat m=S;
    if (S.n_cols != URM.n_cols || S.n_rows != URM.n_rows){
        std::cerr << "ERROR in project_URM : Inconsistent matrices !!!" << std::endl;
    }else{
        for(uword i=0; i<URM.n_rows; ++i){
            for(uword j=0; j<URM.n_cols; ++j){
                if (URM(i,j) != 0){
                    m(i,j)=URM(i,j);
                }
            }
        }
    }
    return m;
}


double evaluate_Obj_Function(const SpMat<double>& URM,const Mat<double>& U,
							 const Mat<double>& H,const SpMat<double>& V,
							 const Mat<double>& U_old,const Mat<double>& H_old,
							 const SpMat<double>& V_old, const double lambda)
{

	mat obj = mat("0");

	for ( std::size_t i(0); i < URM.n_rows ; i++){
		for ( std::size_t j(0); j < URM.n_cols; j++){
			obj = obj +  pow( build_S( i, j, URM, U_old, H_old, V_old) - U.row(i)*H*V.col(j) , 2 );
		}
	}

	return as_scalar( obj + lambda * ( pow( norm (U,"fro"), 2) + pow( norm (H,"fro"), 2) ) );
}




void get_Positive_Matrix(Mat<double> &U)
{
	for (auto i = U.begin() ; i != U.end() ; ++i)
	{
		if (*i < 0)
			*i = 0;
	}
}

uvec get_Vector_Of_Indices(const umat &L){

    // Functino to get de vector of indices from the location matrix
    if (L.n_rows != 2){
        std::cerr << "L has to be a location matrix (size 2xN) !!!" << std::endl;
    }

    uword N_rows=L.row(0).max()+1;
    uvec v(L.n_cols);
    for(uword j=0; j<L.n_cols;++j){
        v(j)=N_rows*L(1,j)+L(0,j);
    }
    return v;
}

///////////////////////////////////
////// DEBUGGING FUNCTIONS ////////
///////////////////////////////////


bool check_V_Constraint(const SpMat<double>& V,const sp_umat &X)
{
	// CHECK 1
	// All zero elements of X must be zero elements of V 

	if (V.n_cols != X.n_cols || V.n_rows != X.n_rows)
	{
		std::cerr << "Inconsistent matrices !!!" << std::endl;
	}

	// I must loop over all the zero elements of X (it's not enough to loop
	// over the non-zero elements of X unfortunately)
	// TODO This loop is inefficient, we could use the CSC format to perform
	// this check, not very important though

    for (std::size_t i = 0 ; i<X.n_rows ; ++i)
	{
        for (std::size_t j = 0 ; j < X.n_cols ; ++j)
		{	
			if (X(i,j)==0 && V(i,j)!=0)
			{
				#ifndef NDEBUG
				std::cerr << "Mismatch found at position (" << i << "," << j << ") " << std::endl;
				#endif
				return false;
			}
		}
	}

	#ifndef NDEBUG
	std::cout << "Consistency with ICM passed!" << std::endl;
	#endif

	// CHECK 2
	// For every column the sum of the elements must be equal to 1
	double eps = 1e-4 ;

	// Loop on each column
	for (std::size_t n=0 ; n< V.n_cols ; n++)
	{
		double sum(0);

		// Loop on the iterator for element access

		for (SpMat<double>::const_iterator i = V.begin_col(n) ; i != V.end_col(n) ; ++i)
		{
			sum += *i ;
		}


		if (sum > 1 + eps || sum < 1 - eps)
		{
			#ifndef NDEBUG
			std::cerr << "Sum is not consistent on column " << n << std::endl;
			#endif
			return false;
		}

	}

	#ifndef NDEBUG
	std::cout << "Consistency with probabilistic interpretation passed! " << std::endl;
	#endif

	return true;
}


bool check_Positive_Matrix(const Mat<double> &U)
{
	// Loop on the elements of the matrix to check positivity 

	for (Mat<double>::const_iterator i = U.begin() ; i != U.end() ; ++i )
	{
		if ( (*i)<0 )
		{
			#ifndef NDEBUG
			std::cerr << "Found negative element!" << std::endl;
			#endif
			return false;
		}
	}

	return true;
}



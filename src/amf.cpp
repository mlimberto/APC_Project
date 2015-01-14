//
//  amf.cpp
//  
//
//
//

#include "amf.h"

//////////////////////////
////  CONSTRUCTORS  //////
//////////////////////////

// Most basic constructor
AMF::AMF() :
lambda_(0),n_max_iter_(100),toll_(0.001)
{
	std::cout << "WARNING : The default constructor has been called, so your instance is still basically empty." << std::endl;
    
}

AMF::AMF(std::string URM_filename, std::string ICM_filename, std::string param_filename)
{
	std::cout << "Creating instance of AMF" << std::endl;

	if (!initialize_Parameters(param_filename))
	{
		std::cout << "WARNING : Initialization of parameters from file didn't work properly" << std::endl;
	}

	if (!import_Sparse_Matrix<double>(URM_filename,URM_,URM_Location_Matrix_,URM_Values_))
	{
		std::cout << "WARNING : Initialization of URM matrix from file didn't work properly" << std::endl;
	}

	if (!import_Sparse_Matrix<unsigned int>(ICM_filename,ICM_,ICM_Location_Matrix_,ICM_Values_))
	{
		std::cout << "WARNING : Initialization of ICM matrix from file didn't work properly" << std::endl;
	}

	// INITIALIZE ALL THE OTHER STUFF

	n_ = URM_.n_rows;
	m_ = ICM_.n_cols;
	k_ = ICM_.n_rows;

	initialize_matrices();


	std::cout << "...done!" << std::endl;

}



//////////////////////////
////  METHODS       //////
//////////////////////////


// TODO La funzione non verifica che tutti i parametri siano stati 
// inizializzati, sarebbe da fare ma non prioritario
bool AMF::initialize_Parameters(std::string filename)
{
	#ifndef NDEBUG
	std::cout << "Initializing parameters of AMF class from file " << filename << std::endl;
	#endif

	// Define the available options through a std::map object
	static std::map<std::string,unsigned int> options_map; 
	options_map["LAMBDA"] = 0;
	options_map["N_MAX_ITER"] = 1;
	options_map["TOLL"] = 2;
	options_map["N_LATENT_FACTORS"] = 3;
	options_map["N_MAX_ITER_GRADIENT"] = 4;
	options_map["TOLL_GRADIENT"] = 5;
	options_map["GRADIENT_STEP"] = 6;



	// Open and read the parameters file 
	std::ifstream param_file(filename);
	std::string option;

	while (param_file >> option)
	{
		#ifndef NDEBUG
		std::cout << "Setting option " << option << std::endl;
		#endif

		// Check if the option is valid
		std::map<std::string,unsigned int>::const_iterator it = options_map.find(option);

		if (it == options_map.end())
		{
			std::cerr << "ERROR : Check consistency of file " << filename << std::endl;
			return false; 
		}

		// Do the mapping and assign the parameter
		unsigned int mapped_option = options_map[option] ;

		switch(mapped_option)
		{
			case 0 : // Next input is value of lambda
			double lambda ; param_file >> lambda ; set_lambda(lambda);
			break;

			case 1:
			unsigned int n_max_iter ; param_file >> n_max_iter ; 
			set_n_max_iter(n_max_iter) ;
			break;

			case 2: 
			double toll ; param_file >> toll ; set_toll(toll);
			break;

			case 3:
			unsigned int n_lf ; param_file >> n_lf ; 
			set_n_latent_factors(n_lf);
			break;

			case 4:
			unsigned int n_max_iter_gradient ; param_file >>n_max_iter_gradient ; 
			set_n_max_iter_gradient(n_max_iter_gradient) ;
			break;

			case 5: 
			double toll_g ; param_file >> toll_g ; set_toll_gradient(toll_g);
			break;

			case 6: 
			double g_step ; param_file >> g_step ; set_gradient_step(g_step);
			break;

			default : 
			std::cout << "Running default option" << endl;
			break;

		}
	}

	return true;
}


void AMF::initialize_matrices(){

	// Initialize U 
    //std::cout<< "Initializing U_old..." << std::endl;

	// U_old_ = 10*randu<mat>(n_,r_);
    //U_old_ = mat(n_,r_,fill::ones);

	// Initialize H
    //std::cout<< "Initializing H_old..." << std::endl;

    //H_old_ = mat(r_,k_,fill::eye);

    // Initialize V
    std::cout <<"Initializing V_old..."<<std::endl;
    uword n_nonzero(0);
    V_old_=sp_mat(ICM_.n_rows,ICM_.n_cols);
    for (uword j = 0 ; j < ICM_.n_cols ; ++j){
        n_nonzero=0;
        for (uword i = 0 ; i<ICM_.n_rows ; ++i){
            if(ICM_(i,j)>0){
                V_old_(i,j)=1;
                n_nonzero+=1;
            }

        }
        if (n_nonzero!=0){
        V_old_.col(j)=V_old_.col(j)/(n_nonzero);
        }else{
            std::cout<<"ERROR : columnn "<<j<<" of ICM is zero everywhere"<<std::endl;
        }
    }
}


void AMF::solve() 
{

	solve_pg_U();

	// Remember to swap the variables

}

void AMF::solve_With_Log() 
{

	solve_pg_U_With_Log();

	// Remember to swap the variables

}


mat AMF::project_URM_by_column(uword j,const mat &S){
    mat m=S;
    if (S.n_rows != URM_.n_rows){
        std::cerr << "ERROR in project_URM_by_column : Inconsistent matrices !!!" << std::endl;
    }else if(S.n_cols != 1){
        std::cerr << "ERROR in project_URM_by_column : You have to pass a column !!!" << std::endl;

    }else{
        uvec q1=get_Vector_Of_Indices(URM_Location_Matrix_);
        uvec q2=find(q1>=URM_.n_rows*j && q1<URM_.n_rows*(j+1));
        m.elem(q1(q2)-URM_.n_rows*j)=URM_Values_(q2);


    }
    return m;
}

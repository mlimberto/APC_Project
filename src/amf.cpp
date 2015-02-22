//
//  amf.cpp
//  
//
//  Authors : Lara Bombardieri, Andrea Brandoli, Matteo Limberto

#include "amf.h"

using namespace arma;

//////////////////////////
////  CONSTRUCTORS  //////
//////////////////////////

// Most basic constructor
AMF::AMF() :
lambda_(0),n_max_iter_(100)
{
	amf_filename_prefix_ = "";
	std::cout << "WARNING : The default constructor has been called, so your instance is still basically empty." << std::endl;
    
}

AMF::AMF(std::string URM_Tr_filename, std::string ICM_filename, std::string param_filename)
{
	std::cout << "Creating instance of AMF ..." << std::endl << std::endl;

	#ifndef NDEBUG
		std::cout << "Debugging printouts are ON" << std::endl;
	#else
		std::cout << "Debugging printouts are OFF" << std::endl;
	#endif

	#ifdef AMFTIME
		std::cout << "Time measurements are ON" << std::endl;
	#else
		std::cout << "Time measurements are OFF" << std::endl;
	#endif	

	if (!initialize_Parameters(param_filename))
	{
		std::cout << "WARNING : Initialization of parameters from file didn't work properly" << std::endl;
	}

	std::cout << "Importing URM " << std::endl;

	if (!import_Sparse_Matrix<double>(URM_Tr_filename,URM_Tr_,URM_Tr_Location_Matrix_,URM_Tr_Values_))
	{
		std::cout << "WARNING : Initialization of URM matrix from file didn't work properly" << std::endl;
	}

	std::cout << "Importing ICM " << std::endl;

	if (!import_Sparse_Matrix<unsigned int>(ICM_filename,ICM_,ICM_Location_Matrix_,ICM_Values_))
	{
		std::cout << "WARNING : Initialization of ICM matrix from file didn't work properly" << std::endl;
	}

	// Initialize the matrices

	n_ = URM_Tr_.n_rows;
	m_ = ICM_.n_cols;
	k_ = ICM_.n_rows;

	initialize_matrices();

	// Set a consistent filename
	amf_filename_prefix_ = "amf_" + std::to_string(r_) + "_" + std::to_string(lambda_) + "_" ;

	#ifndef NDEBUG
	std::cout << "Setting filename name " << amf_filename_prefix_ << std::endl;
	#endif

	// Print some useful information

	std::cout << "N = " << n_ << "  [number of users]" << std::endl;
	std::cout << "M = " << m_ << "  [number of items]" << std::endl;
	std::cout << "L = " << k_ << "  [number of labels]" << std::endl;
	std::cout << "R = " << r_ << "  [number of latent factors]" << std::endl;

}


//////////////////////////
////  METHODS       //////
//////////////////////////


bool AMF::initialize_Parameters(std::string filename)
{
	#ifndef NDEBUG
	std::cout << "Initializing parameters of AMF class from file " << filename << std::endl;
	#endif

	// Define the available options through a std::map object
	static std::map<std::string,unsigned int> options_map; 
	options_map["LAMBDA"] = 0;
	options_map["N_MAX_ITER"] = 1;
	options_map["N_LATENT_FACTORS"] = 3;
	options_map["N_MAX_ITER_GRADIENT"] = 4;
	options_map["TOLL_GRADIENT"] = 5;



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

			default : 
			std::cout << "Running default option" << endl;
			break;

		}
	}

	return true;
}


void AMF::initialize_matrices(){

	#ifndef NDEBUG
	std::cout << "Initializing matrices U,H,V ... " << std::endl;
	#endif

	// Random initialization of matrices

	U_old_ = randu<mat>(n_,r_);

    H_old_ = randu<mat>(r_,k_);

    // Initialize V

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
	for (unsigned int i=0 ; i<n_max_iter_ ; ++i)
	{

		std::cout << std::endl <<"##############" << std::endl << "Iteration " << i+1 << std::endl << "##############" << std::endl << std::endl ;

		std::cout << "SOLVING FOR U ..." << std::endl;

		solve_pg_U();

		std::cout << "SOLVING FOR H ..." << std::endl;

		solve_pg_H();

		std::cout << "SOLVING FOR V ..." << std::endl;

		solve_V();

		std::swap(U_,U_old_);
		std::swap(H_,H_old_);
		std::swap(V_,V_old_);

	}


}

void AMF::solve_With_Log() 
{
	total_logfile_.open(amf_filename_prefix_+"log_iterations.txt");

	for (unsigned int i=0 ; i<n_max_iter_ ; ++i)
	{

		std::cout << std::endl <<"##############" << std::endl << "Iteration " << i+1 << std::endl << "##############" << std::endl << std::endl ;

		std::cout << "SOLVING FOR U ..." << std::endl;

		solve_pg_U_With_Log();

		std::cout << "SOLVING FOR H ..." << std::endl;

		solve_pg_H_With_Log();

		std::cout << "SOLVING FOR V ..." << std::endl;

		solve_V_With_Log();

		std::swap(U_,U_old_);
		std::swap(H_,H_old_);
		std::swap(V_,V_old_);

		total_logfile_ << evaluate_Obj_Function(URM_Tr_,U_old_,H_old_,V_old_,U_old_,H_old_,V_old_,lambda_) << "\n" ;


	}

	total_logfile_.close();

}

void AMF::solve_For_Tuning(std::string mfilename)
{
	// Load validation matrix 
	#ifndef NDEBUG
	std::cout << "Importing URM " << std::endl;
	#endif

	if (!import_Sparse_Matrix<double>(mfilename,URM_Val_,URM_Val_Location_Matrix_,URM_Val_Values_))
	{
		std::cout << "WARNING : Initialization of URM matrix from file didn't work properly" << std::endl;
	}

	#ifndef NDEBUG
	std::cout << "Validation matrix successfully loaded" << std::endl;
	#endif	

	// Run the solver
	std::ofstream val_logfile;
	val_logfile.open(amf_filename_prefix_+"log_validation.txt");
	total_logfile_.open(amf_filename_prefix_+"log_iterations.txt");



	for (unsigned int i=0 ; i<n_max_iter_ ; ++i)
	{

		std::cout << std::endl <<"##############" << std::endl << "Iteration " << i+1 << std::endl << "##############" << std::endl << std::endl ;

		std::cout << "SOLVING FOR U ..." << std::endl;

		solve_pg_U();

		std::cout << "SOLVING FOR H ..." << std::endl;

		solve_pg_H();

		std::cout << "SOLVING FOR V ..." << std::endl;

		solve_V();

		std::swap(U_,U_old_);
		std::swap(H_,H_old_);
		std::swap(V_,V_old_);
		
		double err_val = evaluate_Against_URM_Validation();
		val_logfile << err_val << "\n";

		double obj_fun = evaluate_Obj_Function(URM_Tr_,U_old_,H_old_,V_old_,U_old_,H_old_,V_old_,lambda_) ;
		total_logfile_ << obj_fun << "\n" ;

		std::cout << "Objective function : " << obj_fun << std::endl; 
		std::cout << "RMSE Validation : " << err_val << std::endl;

		// Backup matrices for safety 

		if (i == 0 || i == 10 || i == 20 || i == 30)
		{
			U_old_.save(amf_filename_prefix_+"U"+std::to_string(i));

			H_old_.save(amf_filename_prefix_+"H"+std::to_string(i));

			V_old_.save(amf_filename_prefix_+"V"+std::to_string(i));
		}

	}

	total_logfile_.close();
	val_logfile.close();


}

void AMF::solve_For_Tuning_With_Log(std::string mfilename)
{
	// Load validation matrix 
	#ifndef NDEBUG
	std::cout << "Importing URM " << std::endl;
	#endif

	if (!import_Sparse_Matrix<double>(mfilename,URM_Val_,URM_Val_Location_Matrix_,URM_Val_Values_))
	{
		std::cout << "WARNING : Initialization of URM matrix from file didn't work properly" << std::endl;
	}

	#ifndef NDEBUG
	std::cout << "Validation matrix successfully loaded" << std::endl;
	#endif	

	// Run the solver
	std::ofstream val_logfile;
	val_logfile.open(amf_filename_prefix_+"log_validation.txt");
	total_logfile_.open(amf_filename_prefix_+"log_iterations.txt");



	for (unsigned int i=0 ; i<n_max_iter_ ; ++i)
	{

		std::cout << std::endl <<"##############" << std::endl << "Iteration " << i+1 << std::endl << "##############" << std::endl << std::endl ;

		std::cout << "SOLVING FOR U ..." << std::endl;

		solve_pg_U_With_Log();

		std::cout << "SOLVING FOR H ..." << std::endl;

		solve_pg_H_With_Log();

		std::cout << "SOLVING FOR V ..." << std::endl;

		solve_V_With_Log();

		std::swap(U_,U_old_);
		std::swap(H_,H_old_);
		std::swap(V_,V_old_);

		total_logfile_ << evaluate_Obj_Function(URM_Tr_,U_old_,H_old_,V_old_,U_old_,H_old_,V_old_,lambda_) << "\n" ;
		
		double err_val = evaluate_Against_URM_Validation();
		val_logfile << err_val << "\n";

		#ifndef NDEBUG
		std::cout << err_val << std::endl;
		#endif

	}

	total_logfile_.close();
	val_logfile.close();


}


umat AMF::get_TopN_Recommendation(uword N,bool export_to_file)
{
	umat REC = zeros<umat>(n_,N); // Initialize recommendation matrix

	for (uword i=0 ; i < n_ ; i++)
	{
		// for every row of URM_Tr set a row vector v in which we fill the zeros
		// of URM_Tr with the predictions. If URM_Tr is different from zero,
		// we set v to 0, so we don't consider the evaluations already done
		// from a user

		rowvec v = zeros<rowvec>(m_);
		for (uword j = 0 ; j < m_ ; j++)
		{
			if (URM_Tr_(i,j) == 0 )
			{
				v(j) = as_scalar(U_old_.row(i)*H_old_*V_old_.col(j));
			}
		}

		// now we create the rowvector of indexes sorted in a decreasing order

		urowvec indexes = stable_sort_index(v,"descent");

		// now we keep only the N best recommendations

		indexes.resize(1,N);

		// finally we fill the REC matrix with the N best sorted indexes

		REC.row(i) = indexes;
	}

	if (export_to_file)
	{
		std::ofstream top_N_matrix;
		top_N_matrix.open(amf_filename_prefix_+"top_N.txt") ;

		for (uword i=0 ; i < n_ ; ++i)
		{
			for (uword j = 0 ; j < REC.n_cols ; ++j)
			{
				top_N_matrix << REC(i,j) << " " ;
			}
			top_N_matrix << "\n" ;

		}

		top_N_matrix.close();
	}

	return REC;
}

double AMF::evaluate_Against_URM_Validation()
{
	double result(0);

	for (uword ind = 0 ; ind < URM_Val_Location_Matrix_.n_cols ; ++ind)
	{
		uword i = URM_Val_Location_Matrix_(0,ind);
		uword j = URM_Val_Location_Matrix_(1,ind);

		double val = as_scalar(U_old_.row(i)*H_old_*V_old_.col(j) );
		val = (val - URM_Val_Values_(ind));
		result += val*val;
	}

	result = result / URM_Val_Location_Matrix_.n_cols ;

	return sqrt(result);
}

void AMF::export_Results()
{

	U_old_.save(amf_filename_prefix_+"U");

	H_old_.save(amf_filename_prefix_+"H");

	V_old_.save(amf_filename_prefix_+"V");

}

void AMF::import_Results()
{

	U_old_.load(amf_filename_prefix_+"U");

	H_old_.load(amf_filename_prefix_+"H");

	V_old_.load(amf_filename_prefix_+"V");

}







#include <iostream>
#include <fstream>

#include <armadillo>

#include "my_utils.h"

using namespace arma;
using namespace std;

int main(int argc,char** argv)
{
	// Test projection 

	mat M(2,2,fill::zeros);
	M(0,0) = 0.5; 
	M(1,1) = -0.5;

	M.print("M matrix");

	get_Positive_Matrix(M);

	M.print("M matrix projected");

	// Test gradient algorithm


	mat U(5,3,fill::randu) , A(3,5,fill::eye) ;

	ofstream logtxt;

	logtxt.open("log.txt");

	U.print("U matrix");

	A.print("A matrix");

	cout << norm(U*A,"fro") << endl;

	double step = 0.3;

	mat G(5,3,fill::zeros);

	for (int n = 0 ; n< 30 ; ++n)
	{
		cout << "Iteration " << n+1 << endl ;


		// Compute gradient

		for (uword x =0 ; x< U.n_rows ; ++x)
		{
			for (uword y = 0 ; y < U.n_cols ; ++y)
			{
				double val = 0;

				for (uword k = 0 ; k < U.n_cols ; ++k)
					for(uword j = 0 ; j < U.n_cols ; ++j)
						val += U(x,k)*A(k,j)*A(y,j);

				G(x,y) = 2*val;

			}
		}

		// Update U 
		U = U - step*G ;

		// U.print("U matrix");

		cout << norm(U*A,"fro") << endl;

		// Project U 

		get_Positive_Matrix(U);

		// U.print("U matrix");

		cout << norm(U*A,"fro") << endl << endl;

		logtxt << norm(U*A,"fro") << "\n" ;

	}

	logtxt.close();

	// Introduciamo ora una matrice S 

	mat S(5,5,fill::eye);

	S.print("S matrix");

	// Resettiamo U

	U.randu(5,3);

	U.print("U matrix");

	// Resettiamo G

	G.zeros(5,3);


	// Facciamo partire la risoluzione col gradiente

	ofstream logtxt2;

	logtxt2.open("log2.txt");

	cout << norm(S - U*A,"fro") << endl;


	for (int n = 0 ; n< 30 ; ++n)
	{
		cout << "Iteration " << n+1 << endl ;


		// Compute gradient

		for (uword x =0 ; x< U.n_rows ; ++x)
		{
			for (uword y = 0 ; y < U.n_cols ; ++y)
			{
				double val = 0;

				// Parte quadratica
				for (uword k = 0 ; k < U.n_cols ; ++k)
					for(uword j = 0 ; j < U.n_cols ; ++j)
						val += U(x,k)*A(k,j)*A(y,j);

				// Parte lineare

				double val2 = 0 ;

				for (uword k = 0 ; k < A.n_cols ; ++k)
					val2 += S(x,k)*A(y,k);

				G(x,y) = 2*val - 2*val2 ;

			}
		}

		// Update U 
		U = U - step*G ;

		// U.print("U matrix");

		cout << norm(S - U*A,"fro") << endl;

		// Project U 

		get_Positive_Matrix(U);

		// U.print("U matrix");

		cout << norm(S - U*A,"fro") << endl << endl;

		logtxt2 << norm(S - U*A,"fro") << "\n" ;

		// Per osservare l'andamento aprire il file 
		// log.txt con Matlab e fare un plot

	}

	logtxt2.close();





	return 0;
}

















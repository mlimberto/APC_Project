#include "my_utils.h"


using namespace std;

// ricordarsi che gli indici di riga e colonna in c++ partono da 0 !!

int main(int argc,char **argv)
{

	mat A = mat("1 1 1 1; 0 0 1 1; 1 1 1 1");
	cout << norm (A,"fro") << endl;
	cout << norm (A,"fro")*norm (A,"fro")<< endl;
	cout << pow(norm (A,"fro"),2) << endl;
	cout << A << endl;
	cout << "Numero righe:" << A.n_rows << endl;
	cout << "Numero colonne:" << A.n_cols << endl;
	cout << abs(-1) << endl;
	size_t t(5);
	cout << t << endl;
	return 0;
}

#include <iostream> 

#include "test_lib.h"
#include "lib_brando.h"
#include "testlib_lara.h"

int main(int argc,char** argv)
{
	std:: cout << "Hello world" << std::endl ;
	test_function();
	brando_function();
    lara_function();

	return 0;
}

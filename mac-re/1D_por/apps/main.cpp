#include <iostream>
#include "parser.h"
#include "define.h"
#include "funcs.h"
#include "math.h"
#include "iterate.h"
#include "initial.h"
#include "output.h"
#include <string>
#include "numerics.h"


int main(int argc, const char **argv) {

	options	opt;

	create_default_options(&opt);
	int retCode = parse_options(argc, argv, &opt);
	if (0 != retCode) {
		exit(retCode);
	}

	double *r = new double[(int)NGRID+2];
	double *sigma = new double[(int)NGRID+2];
	double *u = new double[(int)NGRID+2];
	double *w = new double[(int)NGRID+2];

	
	Initialize(r,sigma,u);


	double dt;
	dt = TimeStep(r);
	std::cout << "dt: " << dt/2. << " DT: " << DT << std::endl;		

	if (dt/2. > DT) {
		dt = DT;
	} else {
		dt = dt/2.;
	}

	std::cout << "uj dt: " << dt << std:: endl;	

	getchar();

	PrintData(combine_path("./", "init.dat"),r,sigma,u);

	tIntegrate(r,u,w,dt);

    	std::cout << "\nLefutott..." << std::endl;

	CalculateSigma(r,u,sigma);

	PrintData(combine_path("./", "end.dat"),r,sigma,u);
    
	delete[] r;
	delete[] sigma;
	delete[] u;
	delete[] w;

    	return 0;
}

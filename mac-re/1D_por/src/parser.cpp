#include "parser.h"
#include "iostream"
#include "string.h"
#include "math.h"


void create_default_options(options *opt){

	opt->n				 = 1000;
	opt->m0				 = 1.0;
	opt->r0				 = 0.1;
	opt->r1				 = 5.0;
	opt->sigma0			 = 1.0e-4;
	opt->index			 = -1.0/2.0;
	opt->dt				 = 0.1;
	opt->tStop			 = 5000.0 * 2.0 * M_PI;
	opt->alphaParam		 	 = 1.0e-2;
	opt->h				 = 5.0e-2;
	opt->outDir			 = "";
	opt->filename		 	 = "Distribution.txt";
}




int parse_options(int argc, const char **argv, options *opt){

	int i = 1;

	while (i < argc) {
		std::string p = argv[i];

		// Number of intervals
		if (p == "-n") {
			i++;
			opt->n = atoi(argv[i]);
		}
		else if (p == "-dt") {
			i++;
			opt->dt = atof(argv[i]);
		}
		else if (p == "-tStop") {
			i++;
			opt->tStop = atof(argv[i]);
			opt->tStop *= 2.0 * M_PI;
		}
		else if (p == "-r0") {
			i++;
			opt->r0 = atof(argv[i]);
		}
		else if (p == "-r1") {
			i++;
			opt->r1 = atof(argv[i]);
		}
		else if (p == "-sigma0") {
			i++;
			opt->sigma0 = atof(argv[i]);
		}
		else if (p == "-index") {
			i++;
			opt->index = atof(argv[i]);
		}
		else if (p == "-alpha") {
			i++;
			opt->alphaParam = atof(argv[i]);
		}
		else if (p == "-h") {
			i++;
			opt->h = atof(argv[i]);
		}
		else if (p == "-m0") {
			i++;
			opt->m0 = atof(argv[i]);
		}
		// Print-out location
		else if (p == "-o")	{
			i++;
			opt->outDir = argv[i];
		}
		// Input file
		else if (p == "-f")	{
			i++;
			opt->filename = argv[i];
		}
		else {
			std::cerr << "Invalid switch on command-line: " << p << "." << std::endl;
			return 1;
		}
		i++;
	}
	return 0;
}


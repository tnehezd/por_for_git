//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>

#define ngrid		1000
#define rmin		0.1
#define rmax		5.0
#define dgrid		(rmax-rmin)/(ngrid-1)	//dr

//#define sigma0 0.0001
#define sdexp		-0.5	//surface density profile exponent 
#define alpha_visc	0.01
#define asp_ratio	0.05
#define G			1.0		//ekkor r=1-nél a periódus 2*pi
#define STAR		1.0

using namespace std;

typedef struct options 
{
	//! Number of intervals
	int		n;
	//! Timestep of the integrator
	double	dt;
	//! Timespan of the integration
	double	tStop;
	//! Inner edge of the gaseous disk
	double	r0;
	//! Outer edge of the gaseous disk
	double	r1;
	//! Surface density of the gaseous disk
	double	sigma0;
	//! The index of the profile of the initial surface density of the gaseous disk
	double	index;
	//! The alpha parameter (Shakura & Sunyaev (1973))
	double	alphaParam;
	//! The aspect ratio of the disk (h = H/r)
	double	h;
	//! The mass of the central star
	double	m0;
	//! The output directory where the results will be stored
	string	outDir;
	//! The name of the file where the result will be stored
	string	filename;
} options_t;

string combine_path(string dir, string filename)
{
	if (dir.size() > 0) {
		if (*(dir.end() - 1) != '/' && *(dir.end() - 1) != '\\') {
			return dir + '/' + filename;
		}
		else {
			return dir + filename;
		}
	}
	else {
		return filename;
	}
}

double viscosity(double r)
{
	double nu;
	double alpha_r;
	double a_mod = 0.01;
	double r_dze_i = 1.0;
	double r_dze_o = 20.;
	double Dr_dze_i = 2.0 * r_dze_i * asp_ratio;
	double Dr_dze_o = 2.0 * r_dze_o * asp_ratio;
	double H = asp_ratio * r; 
	double OM = sqrt(G * G * STAR / r / r / r);

	alpha_r = 1.0 - 0.5 * (1.0 - a_mod) * (tanh((r - r_dze_i) / Dr_dze_i) + tanh((r_dze_o - r) / Dr_dze_o));
	//alpha_r = 1.0 + 0.5 * (1.0 - a_mod) * (-erf((r - r_dze_i)/(Dr_dze_i)) + erf((r - r_dze_o)/(Dr_dze_o)));
	//alpha_r = 1.0;
	nu = alpha_visc*asp_ratio*asp_ratio*G*sqrt(r)*alpha_r;

	return nu;
}

double Coeff_1(double r)
{
	double A;
	A = 3*viscosity(r);
	return A;
}

double Coeff_2(double r)
{
	double B;
	B=9*viscosity(r)/(2*r);
	return B;
}

void set_grid(const options *opt, double dx, double *x)
{
	x[0] = opt->r0;
	for (int i = 1; i < opt->n; i++) {
		x[i] = x[i-1] + dx;
	}
	x[opt->n] = opt->r1;
}

void set_initial_profile(const options *opt, const double *x, double *u)
{
	// u a nu*sigma  
	//*u-val mondom meg, hogy tömb, a for-ban már u[i]-val töltöm fel

	// perem:
	u[0] = 0.0;
	for(int i = 1; i < opt->n; i++) {
		double r = x[i];
		u[i] = opt->sigma0 * viscosity(r) * pow(r, opt->index);
	}
	// perem:
	u[opt->n] = 0.0;
}

void	calculate_sigma(int n, const double *x, const double *u, double *sigma)
{
	for(int i = 0; i < n; i++) {
		double r = x[i];
		sigma[i] = u[i] / viscosity(r);
	}
}

void print_data(const string path, int n, const double *x, const double *y)
{
	static char sep = ' ';

	std::ofstream	output;
	output.open(path.c_str(), std::ios_base::trunc);

	for (int i=1; i < n; i++) {
		output << x[i] << sep << y[i] << endl;
	}

	output.flush();
	output.close();
}

void create_default_options(options *opt)
{
	opt->n				 = 1000;
	opt->m0				 = 1.0;
	opt->r0				 = 0.1;
	opt->r1				 = 5.0;
	opt->sigma0			 = 1.0e-4;
	opt->index			 = -1.0/2.0;
	opt->dt				 = 0.1;
	opt->tStop			 = 5000.0 * 2.0 * M_PI;
	opt->alphaParam		 = 1.0e-2;
	opt->h				 = 5.0e-2;
	opt->outDir			 = "";
	opt->filename		 = "Distribution.txt";
}

int parse_options(int argc, const char **argv, options *opt)
{
	int i = 1;

	while (i < argc) {
		string p = argv[i];

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
			cerr << "Invalid switch on command-line: " << p << "." << endl;
			return 1;
		}
		i++;
	}
	return 0;
}


int main(int argc, const char **argv)
{
//	const string baseDir = "C:\\Work\\Projects\\solaris\\src\\SurfaceDensityEvolutionOfAThinDisk\\TestCases";
	const string baseDir = "/home/buba/Dokumentumok/Egyetem/doktori/Programkódok/C++/";
	const string subDir = "Run1";

	options	opt;

	create_default_options(&opt);
	int retCode = parse_options(argc, argv, &opt);
	if (0 != retCode) {
		exit(retCode);
	}

	double	*r = new double[opt.n + 1];
	double	*u = new double[opt.n + 1];
	double	*w = new double[opt.n + 1];
	double	*sigma = new double[opt.n + 1];

	const double dr = (opt.r1 - opt.r0) / opt.n;
	set_grid(&opt, dr, r);
	set_initial_profile(&opt, r, u);
	calculate_sigma(opt.n, r, u, sigma);
 
	string curDir = combine_path(baseDir, subDir);
	print_data(combine_path(curDir, "InitialDistDZE.txt"), opt.n, r, sigma);

	//opt.dt = 0.9 * ((dr*dr) / (2.0 * viscosity(opt.r1)));
	cerr << "dt: " << opt.dt << endl;
	time_t start = time(NULL);
	double t = 0.0;
	do {
		u[0] = u[opt.n] = 0.0;
		w[0] = w[opt.n] = 0.0;
		for (int i = 1; i < opt.n; i++) {
			double du = Coeff_1(r[i])*(u[i+1]-2*u[i]+u[i-1])/(dr*dr) + Coeff_2(r[i])*(u[i+1]-u[i])/dr;
			w[i] = u[i] + opt.dt * du;
		}
		std::swap(u, w);
		// apply boundary conditions

		t += opt.dt;
		cerr << "t: " << t << endl;
	} while (t <= opt.tStop);
	cout << "Total time: " << time(NULL) - start << " s" << endl;

	calculate_sigma(opt.n, r, u, sigma);
	print_data(combine_path(curDir, "FinalDistDZE.txt"), opt.n, r, sigma);

	delete[] r;
	delete[] u;
	delete[] w;
	delete[] sigma;

	return 0;
}

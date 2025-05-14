#ifndef PARSER_H
#define PARSER_H

#include "iostream"
#include "string.h"
#include "math.h"
struct options{
	//! Number of intervals
	int		n;
	//! Timestep of the integrator
	double		dt;
	//! Timespan of the integration
	double		tStop;
	//! Inner edge of the gaseous disk
	double		r0;
	//! Outer edge of the gaseous disk
	double		r1;
	//! Surface density of the gaseous disk
	double		sigma0;
	//! The index of the profile of the initial surface density of the gaseous disk
	double		index;
	//! The alpha parameter (Shakura & Sunyaev (1973))
	double		alphaParam;
	//! The aspect ratio of the disk (h = H/r)
	double		h;
	//! The mass of the central star
	double		m0;
	//! The output directory where the results will be stored
	std::string	outDir;
	//! The name of the file where the result will be stored
	std::string	filename;
};
void create_default_options(options *opt);
int parse_options(int argc, const char **argv, options *opt);

#endif

#include <iostream>
#include "output.h"
#include <fstream>
#include "define.h"
#include <string>
#include "funcs.h"

std::string combine_path(std::string dir, std::string filename) {

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

void PrintData(const std::string path, double *r, double *sigma, double *u) {

	static char sep = ' ';
//	const char path = "proba.txt";
	std::ofstream	output;
	output.open(path.c_str(), std::ios_base::trunc);

	for (int i=1; i < NGRID+1; i++) {
		output << r[i] << sep << sigma[i] << sep << u[i] << sep << Viscosity(r[i]) << std::endl;
	}

	output.flush();
	output.close();

}

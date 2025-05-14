#include "iterate.h"
#include "define.h"
#include "numerics.h"
#include <iostream>
#include "math.h"
#include <iomanip>

//using namespace std;

void tIntegrate(double *r, double *u, double *w, double dt) {
	

  	double time_MODEL_start, time_MODEL_end;
	double time_MODEL = 0.0;
	int InnerLoopCount = 0;
	double PhysicalTime = 0.0;
	int ProfilingCounter = 0;

	do {

		NumDeriv(r,u,w,dt);
		time_MODEL = time_MODEL + dt;
		PhysicalTime = PhysicalTime + dt;
		InnerLoopCount++;

      		std::cout << "\rt: " << std::nouppercase << std::setiosflags(std::ios::scientific) << PhysicalTime / 2.0 / M_PI << " yr"; 
		std::cout << " dt: " << std::nouppercase << std::setiosflags(std::ios::scientific) << dt/2.0/M_PI << " yr";
      		std::cout << " Nsubs: " << InnerLoopCount;
		std::cout << " IT: " << (double) ((time_MODEL)*1000.0/(double)PROFILINGITNUM) << " sec";
		std::cout << std::flush; 
	
	} while (time_MODEL <= TMAX);

}

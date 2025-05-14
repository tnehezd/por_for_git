#include "define.h"
#include "math.h"

const float G		= 1;
const float MSTAR	= 1;

const float NGRID	= 1000;
const float RMIN	= 0.1;
const float RMAX	= 5.0;
const double DD		= (RMAX-RMIN)/(NGRID-1);
	
const double SIGMA_0 	= 0.0001;
const double SDEXP 	= -0.5;
const float ALPHA	= 0.01;
const float ASPRATIO	= 0.05;

const float TMAX 	= 5000.0 * 2.0 * M_PI;
const float DT		= 0.1;
const int PROFILINGITNUM = 10;

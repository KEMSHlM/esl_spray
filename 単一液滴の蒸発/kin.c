#include <stdio.h>
#include <math.h>
#include "path.h"
#include INC_CC_DEF
#include INC_ALM_H
#include INC_3_DEF

void Compute_Chemical_Kinetics(double *omega,double T, double *Y, double *ci, double rog)
{
	
	int sp;

	for(sp=0;sp<spmax;sp++)
	{
		omega[sp] = 0.0;
	}
	
}
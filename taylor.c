#include <stdio.h>
#include <math.h>
#include "taylor.h"

double tcos (double x){
	return 1.0 - x*x*(0.5 - x*x*(1.0/24.0));
}

double tcos_erro (double x){
	return fabs(1*x*x*x*x*x/120.0);
}

double tsqrt (double x){
	return 1.0 + ((x-1)/2.0)*(1-((x-1)/4.0)*(1-(x-1)/2.0));
}

double tsqrt_erro (double x){
	return fabs((15.0/384.0)*sqrt((double)128)*(x-1)*(x-1)*(x-1)*(x-1));
}

int menor (double x, double y){
	return x <= y ? 1 : 0;
}

double fabs (double x){
	return x >= 0.0 ? x : -x;
}
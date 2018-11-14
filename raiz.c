#include "raiz.h"

double pow(double base, double exp)
{
	if(exp < 1e-10 && exp > -1e-10)
		return 1;
	else if(exp < 0)
		return pow(base,exp+1)/base;
	else
		return pow(base,exp-1)*base;
}

double abs(double x)
{
	if(x>0)
		return x;
	else
		return -x;
}

int bissecao (double a, double b, int p, double (*f) (double x), double *r)
{
	double fa = (*f)(a);
	double fb = (*f)(b);
	double c;
	double fc;
	int n = 0;

	while(abs(b-a)/2 > 0.5*pow(10.0,(double)-p))
	{
		c = (a+b)/(2.0);
		fc = (*f)(c);
		if(abs(fc) < 1e-15)
			break;

		if(fc*fa<0.0)
			b = c;
		else
			a = c;

		n++;
	}

	*r = (a+b)/(2.0);
	return n;
}

int IQI (double x0, double x1, double x2, int p, double (*f) (double x), double* r)
{
	double c, Ac, A;
	double x[3], fx[3], ffx[3];
	int i, n = 0;

	x[0] = x0;
	x[1] = x1;
	x[2] = x2;

	for( i = 0 ; i < 3 ; i++ )
	{
		fx[i] = (*f)(x[i]);
		ffx[i] = fx[i]*fx[i];
	}

	while(abs(fx[2]) > 0.5*pow(10.0,(double)-p))
	{

		if(n>100)
			return 0;

		A = ffx[0]*(fx[1]-fx[2]) + ffx[2]*(fx[0]-fx[1]) + ffx[1]*(fx[2]-fx[0]);

		Ac = ffx[0]*(fx[1]*x[2]-fx[2]*x[1]) + ffx[2]*(fx[0]*x[1]-fx[1]*x[0]) + ffx[1]*(fx[2]*x[0]-fx[0]*x[2]);

		c = Ac/A;

		x[0] = x[1];
		x[1] = x[2];
		x[2] = c;

		n++;

		for( i = 0 ; i < 3 ; i++ )
		{
			fx[i] = (*f)(x[i]);
			ffx[i] = (fx[i])*(fx[i]);
		}
		
	}

	(*r) = c;
	return n;
}

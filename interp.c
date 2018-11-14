/*	Implementação de interp.c
	Programado por Guilherme Dantas
	12/9/18
*/

#include "interp.h"

/* Protótitpo das funções encapsuladas */

double NewtonDD (int n, double *xi, double (*f) (double), int i, int j);

void NewtonCoef (int n, double* xi, double (*f) (double), double* bi)
{
   int i;
   for( i = 0 ; i < n ; i++ )
   {
      bi[i] = NewtonDD(n,xi,f,0,i);
   }
}

void Chebyshev (int n, double a, double b, double* xi)
{
   int i;
   for ( i = 0 ; i < n ; i++ )
   {
      // beta = 2*i + 1
      xi[i] = ((b-a)/2)*cos(PI*(2*i+1)/(2*n)) + (a+b)/2;
   }
}

double NewtonAval (int n, double* xi, double* bi, double x)
{
   double sum = bi[n-1];
   int i;

   for( i = n-2 ; i >= 0 ; i-- )
   {
      sum *= (x-xi[i]);
      sum += bi[i];
   }

   return sum;
}

/* Funções encapsuladas */

double NewtonDD (int n, double *xi, double (*f) (double), int i, int j)
{
   if( i == j )
   {
      return (*f)(xi[i]);
   }
   else
   {
      return (NewtonDD(n,xi,f,i+1,j)-NewtonDD(n,xi,f,i,j-1))/(xi[j] - xi[i]);
   }
}
/*	Implementa��o de interp.c
	Programado por Guilherme Dantas
	12/9/18
*/

#include <stdio.h>
#include <math.h>

#define PI (double) 3.1415926535897932384626433832795

/* Calcula as n amostras de Chebyshev para a aproxima��o de uma fun��o qualquer
   dentro do intervalo [a,b], armazenando-as no vetor de double xi */
void Chebyshev (int n, double a, double b, double* xi);

/* Calcula os coeficientes do polin�mio interpolante por diferen�as divididas de Newton,
   de grau n, com amostras xi, da fun��o f, e armazenando-os no vetor de double bi*/
void NewtonCoef (int n, double* xi, double (*f) (double), double* bi);

/* Avalia o polin�mio interpolante por diveren�as divididas de Newton de grau n, com
   amostras xi, coeficientes bi, e retorna a imagem desse polin�mio em x */
double NewtonAval (int n, double* xi, double* bi, double x);
/*
   pendulo.h
   Programado por Guilherme Dantas
   31/10/2018
*/

#include <stdio.h>

/* Constantes do sistema */
#define PI (double) 3.141592653589793238
#define GRAVITY (double) 9.8
#define LENGTH (double) 10.0
#define H (double) 1e-3
#define N_REVOLUCOES 10

double pendulo (double t, double h, double* theta, double* w);

double periodo (double theta_0);

double periodo_simplificado (double theta0);

double absoluto( double x );
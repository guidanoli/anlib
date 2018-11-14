/*
   metiter.h
   Programado por Guilherme Dantas
   07/11/2018
*/

#include <stdio.h>

int Jacobi (int n, double** A, double* b, double* x, double tol);
int GaussSeidel (int n, double** A, double* b, double* x, double tol);
int SOR (int n, double** A, double* b, double* x, double tol, double w);
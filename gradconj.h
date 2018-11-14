/*
   gradconj.h
   Programado por Guilherme Dantas
   14/11/2018
*/

#include <stdio.h>

int GradConj (int n, double** A, double* b, double* x, double tol);
int GradConjJacobi (int n, double** A, double* b, double* x, double tol);
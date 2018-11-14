/*
   qr.h
   Programado por Guilherme Dantas
   26/09/2018 - início da implementação
   27/09/2018 - versão 1.0
*/

#include "matriz.h"

/* Executa a fatoração QR por ortogonalizãção de Gram-Schimdt
   modificado da matriz A m por n, armazenando os valores nas
   matrizes já alocadas Q m por n e R n por n de valores double */
void QR (int m, int n, double** A, double** Q, double** R);

double mmqQR (int m, int n, double** A, double* b, double* x);

void printv( double *v , int n );

void printm( double **M , int m , int n );
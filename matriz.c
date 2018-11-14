/*	Implementação de matriz.h
	Programado por Guilherme Dantas
	15/8/18
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matriz.h"

/* Cria um vetor de dimensão n */
double* vetcria (int n){
	double* new_v = (double *) malloc(sizeof(double)*n);
	return new_v;
}

/* Libera o espaço na memória ocupado por um vetor v */
void vetlibera (double* v){
	free(v);
}

/* Retorna o produto escalar entre dois vetores v e w de dimensão n */
double escalar (int n, double* v, double* w){
	int i;
	double sum = 0;
	for( i = 0 ; i < n ; i ++ )
		sum += v[i]*w[i];
	return sum;
}

/* Retorna a norma-2 de um vetor v de dimensão n */
double norma2 (int n, double* v){
	return sqrt(escalar(n,v,v));
}

/* Retorna uma matriz de dimensões m por n */
double** matcria (int m, int n){
	int i;
	double **colunas = (double **) malloc(sizeof(double *)*m);
	double *linha;
	for( i = 0; i < m; i++ ){
		linha = (double *) malloc(sizeof(double)*n);
		colunas[i] = linha;
	}
	return colunas;
}

/* Libera o espaço na memória ocupado por uma matriz A com m linhas */
void matlibera (int m, double** A){
	int i;
	for( i = 0 ; i < m ; i++ )
		free(A[i]);
	free(A);
}

/* Transpõe uma matriz A m por n e armazena na matriz T n por m */
void transposta (int m, int n, double** A, double** T){
	int i, j;
	for( i = 0 ; i < m ; i++ )
		for( j = 0 ; j < n ; j++ )
			T[j][i] = A[i][j];
}

/* Multiplica uma matriz A m por n por um vetor v de dimensão n e armazena o produto no vetor w */
void multmv (int m, int n, double** A, double* v, double* w){
	int i;
	for( i = 0 ; i < m ; i++ )
		w[i] = escalar(n,A[i],v);
}

/* Multiplica uma matriz A (m por n) por uma matriz B (n por q) e armazena o produto na matriz C (m por q) */
void multmm (int m, int n, int q, double** A, double** B, double** C){
	int i, j;
	double **Bt = matcria(q,n);
	transposta(n,q,B,Bt);
	for( i = 0 ; i < m ; i++ ){
		for( j = 0 ; j < q ; j++ ){
			C[i][j] = escalar(n,A[i],Bt[j]);
		}
	}
	matlibera(q,Bt);
}
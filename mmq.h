/*
   Implementação do módulo mmq.h
   Programado por Guilherme Dantas
   19/09/2018
*/

#include "matriz.h"
#include "sistlinear.h"
#include <math.h>

/* MMQ

   Resolve o sistema A(m,n).x(n) = b(m) pelo Método dos Mínimos Quadrados
   Recebe a matriz A de dimensões m por n, o vetor b de dimensão m
   Preenche o vetor x (já alocado) de dimensão n com a solução do MMQ
*/
void MMQ (int m, int n, double** A, double* b, double* x);

/* RMSE

   Calcula o erro RMSE do resíduo do sistema inconsistente A.x = b, sendo
   a matriz A de dimensões m e n, x um vetor gerado pelo MMQ de dimensão n,
   e b um vetor de dimensão n
*/
double RMSE (int m, int n, double** A, double* b, double* x);

/* periodico

   Ajusta um modelo periódico na forma

   v = c0 + c1 * t + c2 * sin(2*pi*t) + c3 * cos(2*pi*t) + c4 * cos(4*pi*t)

   Recebe n pontos (ti,vi) por dois vetores distintos, ambos de dimensão n.
   Preenche o vetor c (já alocado) de dimensão n com a solução do sistema
   inconsistente por MMQ.
   Retorna o RMSE do resíduo do ajuste.
*/
double periodico (int n, double* t, double* v, double* c);

/* func_periodico
   
   Retorna a imagem da função periódica da forma

   v = c0 + c1 * t + c2 * sin(2*pi*t) + c3 * cos(2*pi*t) + c4 * cos(4*pi*t)

   E recebe os coeficientes num vetor c de dimensão n.
*/
double func_periodico (double *c, double t);
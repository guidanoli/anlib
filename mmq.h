/*
   Implementa��o do m�dulo mmq.h
   Programado por Guilherme Dantas
   19/09/2018
*/

#include "matriz.h"
#include "sistlinear.h"
#include <math.h>

/* MMQ

   Resolve o sistema A(m,n).x(n) = b(m) pelo M�todo dos M�nimos Quadrados
   Recebe a matriz A de dimens�es m por n, o vetor b de dimens�o m
   Preenche o vetor x (j� alocado) de dimens�o n com a solu��o do MMQ
*/
void MMQ (int m, int n, double** A, double* b, double* x);

/* RMSE

   Calcula o erro RMSE do res�duo do sistema inconsistente A.x = b, sendo
   a matriz A de dimens�es m e n, x um vetor gerado pelo MMQ de dimens�o n,
   e b um vetor de dimens�o n
*/
double RMSE (int m, int n, double** A, double* b, double* x);

/* periodico

   Ajusta um modelo peri�dico na forma

   v = c0 + c1 * t + c2 * sin(2*pi*t) + c3 * cos(2*pi*t) + c4 * cos(4*pi*t)

   Recebe n pontos (ti,vi) por dois vetores distintos, ambos de dimens�o n.
   Preenche o vetor c (j� alocado) de dimens�o n com a solu��o do sistema
   inconsistente por MMQ.
   Retorna o RMSE do res�duo do ajuste.
*/
double periodico (int n, double* t, double* v, double* c);

/* func_periodico
   
   Retorna a imagem da fun��o peri�dica da forma

   v = c0 + c1 * t + c2 * sin(2*pi*t) + c3 * cos(2*pi*t) + c4 * cos(4*pi*t)

   E recebe os coeficientes num vetor c de dimens�o n.
*/
double func_periodico (double *c, double t);
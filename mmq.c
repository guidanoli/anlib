/*
   Implementa��o do m�dulo mmq.c
   Programado por Guilherme Dantas
   19/09/2018
*/

#include <stdio.h>
#include <stdlib.h>
#include "mmq.h"

#define PI 3.14159265358979323846

/* Prot�tipo das fun��es encapsuladas pelo m�dulo */

static double func_one ( double x );
static double func_x ( double x );
static double func_sin2pix ( double x );
static double func_cos2pix ( double x );
static double func_cos4pix ( double x );

/* Vari�veis globais a esse m�dulo */

static double (*func_vec[5]) (double x) = 
   { func_one, func_x, func_sin2pix, func_cos2pix, func_cos4pix }; 

/* Fun��es exportadas pelo m�dulo */

void MMQ (int m, int n, double** A, double* b, double* x)
{

   double **At, **M, *b_;
   int *p;

   /* S�o alocados os vetores e matrizes */

   At = matcria(n,m);
   M = matcria(n,n);
   b_ = vetcria(n);
   p = (int *) malloc(sizeof(int)*n);
      
   if( !At || !M || !b_ || !p )
   {
      printf("Faltou memoria!\n");
      return;
   } /* if */

   /* A partir daqui nos certificamos que nenhuma matriz ou vetor � nulo */

   transposta(m,n,A,At);
   multmm(n,m,n,At,A,M);
   multmv(n,m,At,b,b_);

   /* Resolve-se o sistema M.x = b_ por fatora��o LU */

   fatoracao(n,M,p);
   substituicao(n,M,p,b_,x);
   
   /* As matrizes e os vetores s�o liberados da mem�ria */

   matlibera(n,At);
   matlibera(n,M);
   vetlibera(b_);
   free(p);

}

double RMSE (int m, int n, double** A, double* b, double* x)
{
   double *b_;
   double r_norma2 = 0;
   int i;
   
   /* Aloca o vetor b_ que armazenar� o produto matricial A.x */

   b_ = vetcria(m);

   if( !b_ )
   {
      printf("Faltou memoria!\n");
      return -1;
   } /* if */
   
   /* A partir daqui nos certificamos que b_ n�o � nulo */

   /* b_ = A.x */

   multmv(m,n,A,x,b_);
   
   /* Calcula-se a norma2 do res�duo */

   for( i = 0 ; i < m ; i++ )
   {
      r_norma2 += pow(b[i] - b_[i],2);
   } /* for */
   
   /* Divide-se por n e faz-se a raiz quadrada deste valor */

   r_norma2 = sqrt( r_norma2 / m );

   return r_norma2;
}

double periodico (int n, double* t, double* v, double* c)
{
   double **A;
   
   /* Vetor com todas as fun��es para cada coluna da matriz A */

   int i, j;

   /* Matriz A � alocada */

   A = matcria(n,5);

   if( !A )
   {
      printf("Faltou memoria!\n");
      return -1;
   } /* if */

   /* A partir daqui sabe-se que a matriz A n�o � nula */

   /* A matriz A � preenchida com as fun��es de 0 a 5 (coluna)
      no vetor de fun��es para cada valor de ti, variando i
      de 0 a n (linha) */

   for( j = 0 ; j < 5; j++ )
   {
      for( i = 0 ; i < n ; i++ )
      {
         A[i][j] = (*func_vec[j])(t[i]);
      }
   }

   /* � resolvido o sistema inconsistente A.c = v por MMQ,
      armazenando a solu��o em c (j� alocado) */

   MMQ(n,5,A,v,c);

   /* Retorna-se o RMSE do res�duo do ajuste */

   return RMSE(n,5,A,v,c);

}

double func_periodico (double *c, double t)
{
   int i;
   double soma = 0;
   for( i = 0 ; i < 5 ; i++ )
   {
      soma += ((*func_vec[i])(t))*c[i];
   }
   return soma;
}

/* Fun��es encapsuladas pelo m�dulo */

double func_one ( double x )
{
   return 1;
}

double func_x ( double x )
{
   return x;
}

double func_sin2pix ( double x )
{
   return sin(2*PI*x);
}

double func_cos2pix ( double x )
{
   return cos(2*PI*x);
}

double func_cos4pix ( double x )
{
   return cos(4*PI*x);
}
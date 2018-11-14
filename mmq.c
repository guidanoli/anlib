/*
   Implementação do módulo mmq.c
   Programado por Guilherme Dantas
   19/09/2018
*/

#include <stdio.h>
#include <stdlib.h>
#include "mmq.h"

#define PI 3.14159265358979323846

/* Protótipo das funções encapsuladas pelo módulo */

static double func_one ( double x );
static double func_x ( double x );
static double func_sin2pix ( double x );
static double func_cos2pix ( double x );
static double func_cos4pix ( double x );

/* Variáveis globais a esse módulo */

static double (*func_vec[5]) (double x) = 
   { func_one, func_x, func_sin2pix, func_cos2pix, func_cos4pix }; 

/* Funções exportadas pelo módulo */

void MMQ (int m, int n, double** A, double* b, double* x)
{

   double **At, **M, *b_;
   int *p;

   /* São alocados os vetores e matrizes */

   At = matcria(n,m);
   M = matcria(n,n);
   b_ = vetcria(n);
   p = (int *) malloc(sizeof(int)*n);
      
   if( !At || !M || !b_ || !p )
   {
      printf("Faltou memoria!\n");
      return;
   } /* if */

   /* A partir daqui nos certificamos que nenhuma matriz ou vetor é nulo */

   transposta(m,n,A,At);
   multmm(n,m,n,At,A,M);
   multmv(n,m,At,b,b_);

   /* Resolve-se o sistema M.x = b_ por fatoração LU */

   fatoracao(n,M,p);
   substituicao(n,M,p,b_,x);
   
   /* As matrizes e os vetores são liberados da memória */

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
   
   /* Aloca o vetor b_ que armazenará o produto matricial A.x */

   b_ = vetcria(m);

   if( !b_ )
   {
      printf("Faltou memoria!\n");
      return -1;
   } /* if */
   
   /* A partir daqui nos certificamos que b_ não é nulo */

   /* b_ = A.x */

   multmv(m,n,A,x,b_);
   
   /* Calcula-se a norma2 do resíduo */

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
   
   /* Vetor com todas as funções para cada coluna da matriz A */

   int i, j;

   /* Matriz A é alocada */

   A = matcria(n,5);

   if( !A )
   {
      printf("Faltou memoria!\n");
      return -1;
   } /* if */

   /* A partir daqui sabe-se que a matriz A não é nula */

   /* A matriz A é preenchida com as funções de 0 a 5 (coluna)
      no vetor de funções para cada valor de ti, variando i
      de 0 a n (linha) */

   for( j = 0 ; j < 5; j++ )
   {
      for( i = 0 ; i < n ; i++ )
      {
         A[i][j] = (*func_vec[j])(t[i]);
      }
   }

   /* É resolvido o sistema inconsistente A.c = v por MMQ,
      armazenando a solução em c (já alocado) */

   MMQ(n,5,A,v,c);

   /* Retorna-se o RMSE do resíduo do ajuste */

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

/* Funções encapsuladas pelo módulo */

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
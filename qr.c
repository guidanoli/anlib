/*
   qr.c
   Programado por Guilherme Dantas
   26/09/2018 - início da implementação
   27/09/2018 - versão 1.0
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "qr.h"

/* Protótipo das funções encapsuladas pelo módulo */

static void fatoracao (int n, double** a, int* p);
static void substituicao (int n, double** a, int* p, double* b, double* x);
static void MMQ (int m, int n, double** A, double* b, double* x);
static double norma2 (double *v, int n);

/* Código das funções exportadas pelo módulo */

void QR (int m, int n, double** A, double** Q, double** R)
{
   int i, j, k;
   double *W = vetcria( m );

   if( A == NULL || Q == NULL || R == NULL )
   {
      printf("ERRO: Matrizes não alocadas\n");
      return;
   }

   for( i = 0 ; i < n ; i++ )
   {
      for( j = i+1 ; j < n ; j++ )
      {
         R[j][i] = 0;
      }
   }

   for( j = 0 ; j < n ; j++ )
   {
      for( k = 0 ; k < m ; k++ )
      {
         W[k] = A[k][j];
      }

      for( i = 0 ; i < j ; i++ )
      {
         double prod_vet = 0;

         for( k = 0 ; k < m ; k++ )
         {
            prod_vet += Q[k][i] * A[k][i];
         }

         R[i][j] = prod_vet;
         
         for( k = 0 ; k < m ; k++ )
         {
            W[k] -= R[i][j] * Q[k][i];
         }

      }

      R[j][j] = norma2( W , m );

      for( k = 0 ; k < m ; k++ )
      {
         Q[k][j] = W[k] / R[j][j];
      }

   }

   vetlibera(W);

}

double mmqQR (int m, int n, double** A, double* b, double* x)
{

   double **Q, **R, **Qt, *b_;

   Q = matcria( m , n );

   if( Q == NULL )
   {
      printf("ERRO: Faltou memoria ao alocar matriz.\n");
      return -1;
   }

   R = matcria( n , n );

   if( R == NULL )
   {
      printf("ERRO: Faltou memoria ao alocar matriz.\n");
      matlibera(m,Q);
      return -1;
   }

   Qt = matcria( n , m );

   if( Qt == NULL )
   {
      printf("ERRO: Faltou memoria ao alocar matriz.\n");
      matlibera(n,R);
      matlibera(m,Q);
      return -1;
   }

   b_ = vetcria(n);

   if( b == NULL )
   {
      printf("ERRO: Faltou memoria ao alocar matriz.\n");
      matlibera(n,Qt);
      matlibera(n,R);
      matlibera(m,Q);
      return -1;
   }

   QR( m, n, A, Q, R );
   transposta(m,n,Q,Qt);
   multmv(n,m,Qt,b,b_);
   MMQ(n,n,R,b_,x);

   vetlibera( b_ );
   matlibera( n , Qt );
   matlibera( m , Q );
   matlibera( n , R );

   return 0;

}

/* Código das funções encapsuladas pelo módulo */

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
      matlibera(n,At);
      matlibera(n,M);
      vetlibera(b_);
      free(p);
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

void fatoracao (int n, double** a, int* p)
{
   int i, j, k, t;
  for (i=0; i<n; ++i)
    p[i] = i;
  for (j=0; j<n-1; ++j) {
    // find pivô
    int m = j;
    for (i=j+1; i<n; ++i)
      if (fabs(a[i][j]) > fabs(a[m][j]))
        m = i;
    // swap lines: j <-> m
    for (k=0; k<n; ++k) {
      double t = a[j][k];
      a[j][k] = a[m][k];
      a[m][k] = t;
    }
    // register permutation
    t = p[j];
    p[j] = p[m];
    p[m] = t;
    // elimination
    for (i=j+1; i<n; ++i) {
      double f = a[i][j]/a[j][j];
      for (k=j+1; k<n; ++k)
        a[i][k] -= f*a[j][k];
      a[i][j] = f;
    }
  }
}

void substituicao (int n, double** a, int* p, double* b, double* x)
{
   int i, j;
  // forward substitution
  for (i=0; i<n; ++i) {
    double s = 0;
    for (j=0; j<i; ++j) 
      s += a[i][j]*x[j];
    x[i] = b[p[i]] - s;
  }
  // backward substitution
  for (i=n-1; i>=0; --i) {
    double s = 0;
    for (j=i+1; j<n; ++j) 
      s += a[i][j]*x[j];
    x[i] = (x[i] - s) / a[i][i];
  }
}

double norma2 (double *v, int n)
{
   int i;
   double soma = 0;
   for( i = 0 ; i < n ; i++ )
   {
      soma += pow( v[i] , 2 );
   }
   return sqrt(soma);
}

void printv( double *v , int n )
{
   int i;
   printf("[ ");
   for( i = 0 ; i < n ; i++ )
   {
      printf( "%f " , v[i] );
   }
   printf("]\n");
}

void printm( double **M , int m , int n )
{
   int i;
   for( i = 0 ; i < m ; i++ )
   {
      printv( M[i] , n );
   }
}
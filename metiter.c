/*
   metiter.c
   Programado por Guilherme Dantas
   07/11/2018
*/

#include <stdio.h>
#include "metiter.h"
#include "matriz.h"

int Jacobi (int n, double** A, double* b, double* x, double tol)
{
   double **Dinv = NULL, **LU = NULL;
   int i, j, num_iteracoes = 0;
   double erro, *novo_x = NULL, *vet_temp = NULL;

   Dinv = matcria(n,n);

   if( Dinv == NULL )
   {
      printf("Faltou memoria!\n");
      return -1;
   } /* if */

   LU = matcria(n,n);

   if( LU == NULL )
   {
      printf("Faltou memoria!\n");
      matlibera(n,Dinv);
      return -1;
   } /* if */

   novo_x = vetcria(n);

   if( novo_x == NULL )
   {
      printf("Faltou memoria!\n");
      matlibera(n,LU);
      matlibera(n,Dinv);
      return -1;
   } /* if */

   vet_temp = vetcria(n);

   if( vet_temp == NULL )
   {
      printf("Faltou memoria!\n");
      vetlibera(novo_x);
      matlibera(n,LU);
      matlibera(n,Dinv);
      return -1;
   } /* if */

   for( i = 0 ; i < n ; i++ )
   {
      for( j = 0 ; j < n ; j++ )
      {
         if( i == j )
         {
            if( A[i][j] == 0 )
            {
               printf("Diagonal possui um zero");
            } /* if */

            LU[i][j] = 0;
            Dinv[i][j] = (1/A[i][j]);
         } /* if */
         else
         {
            LU[i][j] = A[i][j];
            Dinv[i][j] = 0;
         } /* if */

      } /* for */

   } /* for */

   do {

      for( i = 0 ; i < n ; i++ )
         novo_x[i] = x[i];

      multmv(n,n,LU,novo_x,vet_temp);

      for( i = 0 ; i < n ; i++ )
         vet_temp[i] = -vet_temp[i] + b[i];

      multmv(n,n,Dinv,vet_temp,novo_x);

      for( i = 0 ; i < n ; i++ )
      {
         double temp = x[i];
         x[i] = novo_x[i];
         novo_x[i] -= temp;
      } /* for */

      erro = norma2(n,novo_x);

      num_iteracoes++;

   } while( erro > tol ); /* do-while */

   vetlibera(vet_temp);
   vetlibera(novo_x);
   matlibera(n,LU);
   matlibera(n,Dinv);

   return num_iteracoes;

}

int GaussSeidel (int n, double** A, double* b, double* x, double tol)
{
   int i, j, num_iteracoes = 0;
   double erro, erro2, *x_ant = NULL;

   x_ant = vetcria(n);

   if( x_ant == NULL )
   {
      printf("Faltou memoria!\n");
      return -1;
   } /* if */

   do {

      for( i = 0 ; i < n ; i++ )
      {
         double sum = 0;
         
         x_ant[i] = x[i];

         for( j = 0 ; j < n ; j++ )
         {

            if( i != j )
            {
               sum += A[i][j] * x[j];
            } /* if */

         } /* for */

         x[i] = (1/(A[i][i])) * ( b[i] - sum );

         x_ant[i] -= x[i];

      } /* for */

      erro = norma2(n,x_ant);

      num_iteracoes++;

   } while( erro > tol ); /* do-while */

   vetlibera(x_ant);

   return num_iteracoes;
}

int SOR (int n, double** A, double* b, double* x, double tol, double w)
{
   int i, j, num_iteracoes = 0;
   double erro, erro2, *x_ant = NULL, *x_diff = NULL;

   x_ant = vetcria(n);

   if( x_ant == NULL )
   {
      printf("Faltou memoria!\n");
      return -1;
   } /* if */

   x_diff = vetcria(n);

   if( x_diff == NULL )
   {
      printf("Faltou memoria!\n");
      vetlibera(x_ant);
      return -1;
   } /* if */

   do {

      for( i = 0 ; i < n ; i++ )
      {
         double sum = 0;
         
         x_ant[i] = x[i];
         x_diff[i] = x[i];

         for( j = 0 ; j < n ; j++ )
         {

            if( i != j )
            {
               sum += A[i][j] * x[j];
            } /* if */

         } /* for */

         x[i] = (1/(A[i][i])) * ( b[i] - sum );

         x_diff[i] -= x[i];

      } /* for */

      for( i = 0 ; i < n ; i++ )
      {
         x[i] = (1-w)*x_ant[i] + w*x[i];
      } /* for */

      erro = norma2(n,x_diff);

      num_iteracoes++;

   } while( erro > tol ); /* do-while */

   vetlibera(x_ant);

   return num_iteracoes;
}
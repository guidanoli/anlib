/*
   gradconj.c
   Programado por Guilherme Dantas
   14/11/2018
*/

#include "gradconj.h"
#include "matriz.h"

int GradConj (int n, double** A, double* b, double* x, double tol)
{
   double *d = NULL, *r = NULL, *r_ant = NULL, *v_aux = NULL;
   double alpha, beta;

   int i, j, k;
   int num_iteracoes = 0;

   d = vetcria(n);

   if( d == NULL )
   {
      printf("Faltou memoria!");
      return -1;
   } /* if */

   r = vetcria(n);

   if( r == NULL )
   {
      vetlibera(d);
      printf("Faltou memoria!");
      return -1;
   } /* if */

   r_ant = vetcria(n);

   if( r_ant == NULL )
   {
      vetlibera(r);
      vetlibera(d);
      printf("Faltou memoria!");
      return -1;
   } /* if */

   v_aux = vetcria(n);

   if( v_aux == NULL )
   {
      vetlibera(r_ant);
      vetlibera(r);
      vetlibera(d);
      printf("Faltou memoria!");
      return -1;
   } /* if */

   for( i = 0 ; i < n ; i++ )
   {
      double Axi = 0;

      for( j = 0 ; j < n ; j++ )
      {
         Axi += A[i][j] * x[j];
      } /* for */

      d[i] = r_ant[i] = b[i] - Axi;

   } /* for */

   for( k = 0 ; k < n ; k++ )
   {

      double d_aux = 0;

      if( norma2(n,r_ant) < tol )
      {
         break;
      } /* if */

      alpha = 0;
      
      for( i = 0 ; i < n ; i++ )
      {
         alpha += r_ant[i]*r_ant[i];
      } /* for */

      multmv(n,n,A,d,v_aux);

      for( i = 0 ; i < n ; i++ )
      {
         d_aux += v_aux[i] * d[i];
      } /* for */

      alpha /= d_aux;

      for( i = 0 ; i < n ; i++ )
      {
         x[i] += alpha*d[i];

         d_aux = 0;

         for( j = 0 ; j < n ; j++ )
         {
            d_aux += A[i][j]*d[j];
         } /* for */

         r[i] = r_ant[i] - alpha*d_aux;

      } /* for */

      beta = 0;
      d_aux = 0;

      for( i = 0 ; i < n ; i++ )
      {
         beta += r[i]*r[i];
         d_aux += r_ant[i]*r_ant[i];
      } /* for */

      beta /= d_aux;

      for( i = 0 ; i < n ; i++ )
      {
         d[i] = r[i] + beta*d[i];
         r_ant[i] = r[i];
      } /* for */

      num_iteracoes++;
   } /* for */

   vetlibera(r);
   vetlibera(d);
   vetlibera(r_ant);
   vetlibera(v_aux);

   return num_iteracoes;

}

int GradConjJacobi (int n, double** A, double* b, double* x, double tol)
{
   double *d = NULL, *r = NULL, *r_ant = NULL, *v_aux = NULL, *z = NULL, *z_ant = NULL;
   double alpha, beta;

   int i, j, k;
   int num_iteracoes = 0;

   d = vetcria(n);

   if( d == NULL )
   {
      printf("Faltou memoria!");
      return -1;
   } /* if */

   r = vetcria(n);

   if( r == NULL )
   {
      vetlibera(d);
      printf("Faltou memoria!");
      return -1;
   } /* if */

   r_ant = vetcria(n);

   if( r_ant == NULL )
   {
      vetlibera(r);
      vetlibera(d);
      printf("Faltou memoria!");
      return -1;
   } /* if */

   z = vetcria(n);

   if( z == NULL )
   {
      vetlibera(r_ant);
      vetlibera(r);
      vetlibera(d);
      printf("Faltou memoria!");
      return -1;
   } /* if */

   z_ant = vetcria(n);

   if( z_ant == NULL )
   {
      vetlibera(z);
      vetlibera(r_ant);
      vetlibera(r);
      vetlibera(d);
      printf("Faltou memoria!");
      return -1;
   } /* if */

   v_aux = vetcria(n);

   if( v_aux == NULL )
   {
      vetlibera(z_ant);
      vetlibera(z);
      vetlibera(r_ant);
      vetlibera(r);
      vetlibera(d);
      printf("Faltou memoria!");
      return -1;
   } /* if */

   for( i = 0 ; i < n ; i++ )
   {
      double Axi = 0;

      for( j = 0 ; j < n ; j++ )
      {
         Axi += A[i][j] * x[j];
      } /* for */

      r_ant[i] = b[i] - Axi;
      
   } /* for */

   for( i = 0 ; i < n ; i++ )
   {
      d[i] = z_ant[i] = r_ant[i]/A[i][i];

   } /* for */

   for( k = 0 ; k < n ; k++ )
   {

      double d_aux = 0;

      if( norma2(n,r_ant) < tol )
      {
         break;
      } /* if */

      alpha = 0;
      
      for( i = 0 ; i < n ; i++ )
      {
         alpha += r_ant[i]*z_ant[i];
      } /* for */

      multmv(n,n,A,d,v_aux);

      for( i = 0 ; i < n ; i++ )
      {
         d_aux += v_aux[i] * d[i];
      } /* for */

      alpha /= d_aux;

      for( i = 0 ; i < n ; i++ )
      {
         x[i] += alpha*d[i];

         d_aux = 0;

         for( j = 0 ; j < n ; j++ )
         {
            d_aux += A[i][j]*d[j];
         } /* for */

         r[i] = r_ant[i] - alpha*d_aux;
         z[i] = r[i]/A[i][i];

      } /* for */

      beta = 0;
      d_aux = 0;

      for( i = 0 ; i < n ; i++ )
      {
         beta += r[i]*z[i];
         d_aux += r_ant[i]*z_ant[i];
      } /* for */

      beta /= d_aux;

      for( i = 0 ; i < n ; i++ )
      {
         d[i] = z[i] + beta*d[i];
         r_ant[i] = r[i];
         z_ant[i] = z[i];
      } /* for */

      num_iteracoes++;
   } /* for */

   vetlibera(r);
   vetlibera(d);
   vetlibera(r_ant);
   vetlibera(v_aux);
   vetlibera(z_ant);
   vetlibera(z);

   return num_iteracoes;
}
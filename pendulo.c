/*
   pendulo.c
   Programado por Guilherme Dantas
   31/10/2018
*/

#include <stdio.h>
#include <math.h>
#include "pendulo.h"

double absoluto( double x )
{
   if( x > 0 )
   {
      return x;
   }
   else
   {
      return -x;
   }
}

double pendulo (double t, double h, double* theta, double* w)
{
   /* theta" = -(g/l)*sin(theta)
      theta' = w 
      w' = theta" */

   double k_theta[4], k_w[4];
   double new_w, new_theta;
   int i;

   k_theta[0] = h*(*w);
   k_w[0] = -(GRAVITY/LENGTH)*sin((*theta))*h;

   for( i = 1 ; i < 4 ; i++ )
   {
      k_theta[i] = h*( *w + k_w[i-1]/2 );
      k_w[i] = -(GRAVITY/LENGTH)*sin( (*theta) + k_theta[i-1] )*h;
   } /* for */

   new_w = (*w) + (k_w[0] + 2*k_w[1] + 2*k_w[2] + k_w[3])/6;
   new_theta = (*theta) + (k_theta[0] + 2*k_theta[1] + 2*k_theta[2] + k_theta[3])/6;

   *w = new_w;
   *theta = new_theta;

   return t+h;
}

double periodo (double theta_0)
{
   double theta = theta_0, w = 0, t = 0;
   int inversoes = 0;

   double t1, t2, w1, w2;

   while( inversoes < 2*N_REVOLUCOES )
   {
      w1 = w;
      t1 = t;

      t = pendulo( t , H , &theta, &w );

      if( w1*w < 0 || ( w1*w == 0 && ( w1 > w ) ) )
      {
         /* Quando w inverte de sinal ou chega ao ponto de inflexão (raro) */
         inversoes++;
      } /* if */

      w2 = w;
      t2 = t;

   } /* while */

   return (2*(t1+(absoluto(w1)/(absoluto(w1)+absoluto(w2))*(t2-t1))))/(2*N_REVOLUCOES);

}

double periodo_simplificado (double theta0)
{
   return 2*PI*sqrt(LENGTH/GRAVITY);
}
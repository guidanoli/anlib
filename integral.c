/* Implementação do módulo integral
   Programado por Guilherme Dantas
   10/10/2018 */

   #include <stdio.h>
   #include "integral.h"
   #include <math.h>

/* Protótipos das funções encapsuladas */

   static double fabs( double x );

/* Funções exportadas pelo módulo */

   //Derivada

   double derivada (double (*f) (double x), double x, double h)
   {
      if( f == NULL || h == 0 )
      {
         return 0.0;
      }

      return ((*f)(x+h)-(*f)(x-h))/(2*h);
   }

   //H ótimo

   double h_otimo (double (*f) (double x), double (*fl) (double x), double x)
   {
      double h, h_min, e_min;
      int i;

      if( f == NULL || fl == NULL )
      {
         return 0.0;
      }

      h_min = 1e-1;
      e_min = fabs( derivada(f,x,h_min) - (*fl)(x) );

      for( i = 0, h = 1e-1 ; i < 12 ; i++, h /= 10 )
      {
         double e = fabs(derivada(f,x,h) - (*fl)(x));
         if( e < e_min )
         {
            e_min = e;
            h_min = h;
         }
      }

      return h_min;

   }

   //Integração numérica por Simpson

   double simpson (double (*f) (double), double a, double b, int n)
   {

      double f_esq, f_dir, h = (b-a)/n, soma = 0;
      int i;

      if( f == NULL )
      {
         return 0.0;
      }

      f_dir = (*f)(a);

      for( i = 0 ; i < n ; i++ )
      {
         f_esq = f_dir;
         f_dir = (*f)(a+h*(i+1));
         soma += (h/6)*( f_esq + 4*((*f)(a+h*i+(h/2))) + f_dir );
      }

      return soma;
   }

   //Integração númerica por ponto médio

   double pontomedio (double (*f) (double), double a, double b, int n)
   {

      double h = (b-a)/n, soma = 0;
      int i;

      if( f == NULL )
      {
         return 0.0;
      }

      for( i = 0 ; i < n ; i++ )
      {
         soma += h*(*f)(a+i*h+h/2);
      }

      return soma;

   }

/* Funções encapsuladas pelo módulo */

   double fabs( double x )
   {
      if( x < 0 )
         return -x;
      else
         return x;
   }
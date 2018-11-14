/*
   simpson.c
   Programado por Guilherme Dantas
   17/10/2018
*/

#include "simpson.h"
#include <stdio.h>
#include <math.h>

/* Protótipo das funções encapsuladas pelo módulo */

static double Simpson( double a ,
                double b ,
                double (*f) (double x) );

static double TransformaQuadratura( double a ,
                             double b ,
                             double t ,
                             double (*f) (double x) );

/* Funções exportadas pelo módulo */

/*
------------------------------------------
ERRO DE REGRA DE SIMPSON COM DOIS PASSOS
------------------------------------------
*/

double DoubleSimpson ( double a ,
                       double b ,
                       double (*f) (double x) ,
                       double* v )
{
   double m = (a+b)/2;
   double Sab, Sam, Smb;
   double erro;
   Sab = Simpson(a,b,f);
   Sam = Simpson(a,m,f);
   Smb = Simpson(m,b,f);
   erro = (Sab - ( Sam + Smb ))/15;
   (*v) = Sam + Smb - erro;
   return erro;
}

/*
------------------------------------------
REGRA DE SIMPSON ADAPTATIVO RECURSIVA
------------------------------------------
*/

double AdaptiveSimpson ( double a ,
                         double b ,
                         double (*f) (double x) ,
                         double tol )
{
   double m = (a+b)/2;
   double v;
   double delta = DoubleSimpson(a,b,f,&v);
   if( delta > 15*tol )
      return AdaptiveSimpson(a,m,f,tol/2) + AdaptiveSimpson(m,b,f,tol/2);
   else
      return v;
}

/*
------------------------------------------
QUADRATURA COM 2 AMOSTRAS
------------------------------------------
*/

double Quadratura2 ( double a ,
                     double b ,
                     double (*f) (double x) )
{
   double x[2] = { - sqrt((double)1/3) , 
                     sqrt((double)1/3) };

   double c[2] = { 1 , 1 };
   double soma = 0;

   int i;

   for( i = 0 ; i < 2 ; i++ )
      soma += c[i]*TransformaQuadratura(a,b,x[i],f);

   return soma;
}

/*
------------------------------------------
QUADRATURA COM 3 AMOSTRAS
------------------------------------------
*/

double Quadratura3 ( double a ,
                     double b ,
                     double (*f) (double x) )
{
   double x[3] = { - sqrt((double)3/5) , 
                     0 ,
                     sqrt((double)3/5) };

   double c[3] = { (double)5/9 , (double)8/9 , (double)5/9 };
   double soma = 0;

   int i;

   for( i = 0 ; i < 3 ; i++ )
      soma += c[i]*TransformaQuadratura(a,b,x[i],f);

   return soma;
}

/* Corpo das funções encapsuladas pelo módulo */

double Simpson( double a ,
                double b ,
                double (*f) (double x) )
{
   double m = (a+b)/2;
   double h = b-a;
   return (h/6)*((*f)(a)+4*(*f)(m)+(*f)(b));
}

double TransformaQuadratura( double a ,
                             double b ,
                             double t ,
                             double (*f) (double x) )
{
   return (*f)(( (b-a)*t + (b+a) )/2) * (b-a)/2;
}
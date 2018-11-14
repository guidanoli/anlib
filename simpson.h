/*
   simpson.h
   Programado por Guilherme Dantas
   17/10/2018
*/

#include <stdio.h>

/* Calcula a diferen�a da aplica��o da Regra de Simpson para um passo e dois
   passo sobre um mesmo intervalo, retornando 1/15 desta diferen�a, e armazenando
   na vari�vel *v a aplica��o da Regra de Simpson com dois passos mais a diferen�a */
double DoubleSimpson ( double a ,
                       double b ,
                       double (*f) (double x) ,
                       double* v ) ;

/* Calcula e retorna a integral sobre o dom�nio [a,b] da fun��o f, com a toler�ncia
   m�xima "tol", atrav�s da Regra de Simspon Adaptativa */
double AdaptiveSimpson ( double a ,
                         double b ,
                         double (*f) (double x) ,
                         double tol ) ;

/* Calcula e retorna a integral sobre o dom�nio [a,b] da fun��o f atrav�s da
   Quadratura de Gauss com 2 amostras */
double Quadratura2 ( double a ,
                     double b ,
                     double (*f) (double x) ) ;

/* Calcula e retorna a integral sobre o dom�nio [a,b] da fun��o f atrav�s da
   Quadratura de Gauss com 3 amostras */
double Quadratura3 ( double a ,
                     double b ,
                     double (*f) (double x) ) ;
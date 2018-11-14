/*
   simpson.h
   Programado por Guilherme Dantas
   17/10/2018
*/

#include <stdio.h>

/* Calcula a diferença da aplicação da Regra de Simpson para um passo e dois
   passo sobre um mesmo intervalo, retornando 1/15 desta diferença, e armazenando
   na variável *v a aplicação da Regra de Simpson com dois passos mais a diferença */
double DoubleSimpson ( double a ,
                       double b ,
                       double (*f) (double x) ,
                       double* v ) ;

/* Calcula e retorna a integral sobre o domínio [a,b] da função f, com a tolerância
   máxima "tol", através da Regra de Simspon Adaptativa */
double AdaptiveSimpson ( double a ,
                         double b ,
                         double (*f) (double x) ,
                         double tol ) ;

/* Calcula e retorna a integral sobre o domínio [a,b] da função f através da
   Quadratura de Gauss com 2 amostras */
double Quadratura2 ( double a ,
                     double b ,
                     double (*f) (double x) ) ;

/* Calcula e retorna a integral sobre o domínio [a,b] da função f através da
   Quadratura de Gauss com 3 amostras */
double Quadratura3 ( double a ,
                     double b ,
                     double (*f) (double x) ) ;
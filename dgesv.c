#include "dgesv.h"
#include <stdio.h>
#include <math.h>

int my_dgesv(int n, int nrhs, double *a, double *b)
{
  /* add/change the arguments according to your implementation needs */

  // [...] Write your implementation here

  // Definicion de macros para notacion de matrices
  #define a(i,j) a[(i)*n + (j)]
  #define b(i,j) b[(i)*nrhs + (j)]


  // Verificar si tiene sentido el sistema de ecuaciones y reportar error
  if (n <= 0 || nrhs <= 0 || a == NULL || b == NULL)
  return -1;

  // Main: Eliminacion hacia adelante con pivoteo parcial

    for (int k = 0; k < n - 1; k++) {
        // SubMain_1: Buscar fila con maximo valor absoluto en la columna k para evitar 0 en el pivote
          int pivote_fila = k;
          double max_val = fabs(a(k,k));
            for (int i = k + 1; i < n; i++) {
              if (fabs(a(i,k)) > max_val) {
                max_val = fabs(a(i,k));
                pivote_fila = i;
              }
            }

        // SubMain_2: Verificar error si el pivote es muy cercano a cero
            double tol = 1e-8;
            if (fabs(a(pivote_fila,k)) < tol){
              printf("Error! El pivote de la columna %d está muy cerca de 0\n", k);
              return -1;
            }

        // SubMain_3: Intercambiamos filas en la matriz del sistema si el pivote no esta en la fila k
            if (pivote_fila != k){
              for (int j=k;j<n;j++){
                double tmp = a(k,j);
                a(k,j)=a(pivote_fila,j);
                a(pivote_fila,j)=tmp;
              }
            //Lo mismo para la matriz b usando nrhs=number of rigth hand size 
              for (int j=0;j<nrhs;j++){
                double tmp = b(k,j);
                b(k,j)=b(pivote_fila,j);
                b(pivote_fila,j)=tmp;
              }
            }
         
        // SubMain_4: Algoritmo de eliminacion de Gauss
        for (int i = k + 1; i < n; i++){
          double mik = a(i,k) / a(k,k);
            a(i,k) = 0.0;
            for (int j = k+1; j<n;j++){
              a(i,j) = a(i,j)-mik * a(k,j);
            } 
            for (int j = 0; j<nrhs;j++){
              b(i,j) = b(i,j)-mik * b(k,j);
            } 
        }
        
        // Substitucion hacia atras
        for (int col = 0; col < nrhs; col++) {   // bucle del vector b (normal nrhs=1)
          for (int i = n - 1; i >= 0; i--) {     // Ultima fila del sistema
            double sum = b(i,col);              // Guardamos en sum el termino ind de b 
              for (int j = i + 1; j < n; j++) {
                sum = sum - a(i,j) * b(j,col);   // Sustituye (cuando entra el bucle por el anterior)
              }
            b(i,col) = sum / a(i,i); // divido por su coeficiente de la matriz 
          }
        }
        
      }

  return 0;
}

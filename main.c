#include "matrixAlgorithms.h"
#include <stdio.h>
int main(){
    //testing
    for(int n = 10; n < 100000; n*=10 ){
        printf("n = %d\n",n);
        matrix* A = hilbertMatrix(n);
        matrix** QRc = classicalGramSchmidt(A);
        matrix** QRm = modifiedGramSchmidt(A);
        matrix* temp = matrixMultiply(transpose(QRc[0]), QRc[0]);
        freeMatrix(QRc[0]);
        freeMatrix(QRc[1]);
        double cgsRemainder = frobeniusNorm(matrixSubtract(identityMatrix(n,n),temp));
        temp = matrixMultiply(transpose(QRm[0]), QRm[0]);
        double mgsRemainder = frobeniusNorm(matrixSubtract(identityMatrix(n,n),temp));
        printf("CGS Remainder: %f\n",cgsRemainder);
        printf("MGS Remainder: %f\n",mgsRemainder);
    }
}
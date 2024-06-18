#include "vector.h"
#pragma once
typedef struct
{
    int rows;
    int cols;
    double** data;
} matrix;

matrix* makeMatrix(int rows, int cols);
matrix* copyMatrix(matrix* A);
matrix* matrixMultiply(matrix* A, matrix* B);
void printMatrix(matrix* A);
void freeMatrix(matrix* A);
void fillMatrixColumn(matrix* A, int col, double* data);
matrix* transpose(matrix* A);
matrix* randomMatrix(int rows, int cols);
matrix* randomWholeNumberMatrix(int rows, int cols, int max);
matrix* identityMatrix(int rows, int cols);
matrix* matrixAdd(matrix* A, matrix* B);
matrix* matrixSubtract(matrix* A, matrix* B);
matrix* scalarMultiplyMatrix(matrix* A, double scalar);
vector *getRow(matrix *A, int row);
vector *getColumn(matrix *A, int col);
int matrixEquals(matrix *A, matrix *B);
matrix* hilbertMatrix(int n);
double frobeniusNorm(matrix* A);
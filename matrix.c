//This file contains the code for functionality of a matrix like object
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
matrix *makeMatrix(int rows, int cols)
{
    matrix *m = malloc(sizeof(matrix));
    m->rows = rows;
    m->cols = cols;
    m->data = (double **) malloc(sizeof(double *) * rows);
    for(int i = 0; i < rows; i++){
        m->data[i] = (double *) malloc(sizeof(double) * cols);
    }
    return m;
}

matrix *copyMatrix(matrix *A)
{
    matrix *B = makeMatrix(A->rows, A->cols);
    for (int i = 0; i < A->rows; i++)
    {
        for (int j = 0; j < A->cols; j++)
        {
            B->data[i][j] = A->data[i][j];
        }
    }
    return B;
}
matrix *matrixMultiply(matrix *A, matrix *B)
{
    if (A->cols != B->rows)
    {
        printf("Error: Matrix dimensions do not match\n");
        return NULL;
    }
    matrix *C = makeMatrix(A->rows, B->cols);
    for (int i = 0; i < A->rows; i++)
    {
        for (int j = 0; j < B->cols; j++)
        {
            double sum = 0;
            for (int k = 0; k < B->rows; k++)
            {
                sum += A->data[i][k] * B->data[k][j];
            }
            C->data[i][j] = sum;
        }
    }
    return C;
}

void printMatrix(matrix *A)
{
    for (int i = 0; i < A->rows; i++)
    {
        printf("| ");
        for (int j = 0; j < A->cols; j++)
        {   
            printf("%f ", A->data[i][j]);
        }
        printf("|");
        printf("\n");
    }
}

void freeMatrix(matrix *A)
{
    for (int i = 0; i < A->rows; i++)
    {
        free(A->data[i]);
    }
    free(A->data);
    free(A);
}

void fillMatrixColumn(matrix *A, int col, double *data)
{
    for (int i = 0; i < A->rows; i++)
    {
        A->data[i][col] = data[i];
    }
}

matrix *transpose(matrix *A)
{
    matrix *B = makeMatrix(A->cols, A->rows);
    for (int i = 0; i < A->rows; i++)
    {
        for (int j = 0; j < A->cols; j++)
        {
            B->data[j][i] = A->data[i][j];
        }
    }
    return B;
}
vector *getColumn(matrix *A, int col)
{
    double *data = malloc(sizeof(double) * A->rows);
    for (int i = 0; i < A->rows; i++)
    {
        data[i] = A->data[i][col];
    }
    return makeVectorFromData(A->rows, data);
}
vector *getRow(matrix *A, int row)
{
    double *data = malloc(sizeof(double) * A->cols);
    for (int i = 0; i < A->cols; i++)
    {
        data[i] = A->data[row][i];
    }
    return makeVectorFromData(A->cols, data);
}
matrix* randomMatrix(int rows, int cols){
    srand(time(NULL));
    matrix* A = makeMatrix(rows, cols);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            A->data[i][j] = (double)rand() / (double)RAND_MAX;
        }
    }
    return A;
}
matrix* randomWholeNumberMatrix(int rows, int cols, int max){
    srand(time(NULL));
    matrix* A = makeMatrix(rows, cols);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            A->data[i][j] = rand() % max;
        }
    }
    return A;
}
matrix* identityMatrix(int rows, int cols){
    matrix* A = makeMatrix(rows, cols);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            A->data[i][j] = (i == j) ? 1 : 0;
        }
    }
    return A;
}

matrix* matrixAdd(matrix* A, matrix* B){
    if(A->rows != B->rows || A->cols != B->cols){
        printf("Error: Matrix dimensions do not match\n");
        return NULL;
    }
    matrix* C = makeMatrix(A->rows, A->cols);
    for(int i = 0; i < A->rows; i++){
        for(int j = 0; j < A->cols; j++){
            C->data[i][j] = A->data[i][j] + B->data[i][j];
        }
    }
    return C;
}

matrix* matrixSubtract(matrix* A, matrix* B){
    if(A->rows != B->rows || A->cols != B->cols){
        printf("Error: Matrix dimensions do not match\n");
        return NULL;
    }
    matrix* C = makeMatrix(A->rows, A->cols);
    for(int i = 0; i < A->rows; i++){
        for(int j = 0; j < A->cols; j++){
            C->data[i][j] = A->data[i][j] - B->data[i][j];
        }
    }
    return C;
}

matrix* scalarMultiplyMatrix(matrix* A, double scalar){
    matrix* B = makeMatrix(A->rows, A->cols);
    for(int i = 0; i < A->rows; i++){
        for(int j = 0; j < A->cols; j++){
            B->data[i][j] = A->data[i][j] * scalar;
        }
    }
    return B;
}

int matrixEquals(matrix* A, matrix* B){
    if(A->rows != B->rows || A->cols != B->cols){
        return 0;
    }
    for(int i = 0; i < A->rows; i++){
        for(int j = 0; j < A->cols; j++){
            if(A->data[i][j] != B->data[i][j]){
                return 0;
            }
        }
    }
    return 1;
}

matrix* hilbertMatrix(int n){
    matrix* A = makeMatrix(n, n);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            A->data[i][j] = 1.0 / (i + j + 1);
        }
    }
    return A;
}

double frobeniusNorm(matrix* A){
    double sum = 0;
    for(int i = 0; i < A->rows; i++){
        for(int j = 0; j < A->cols; j++){
            sum += (A->data[i][j])*(A->data[i][j]);
        }
    }
    return sqrt(sum);
}
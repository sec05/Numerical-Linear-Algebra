// This file contains algorithms to factorize matricies or solve linear systems
#include "matrixAlgorithms.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

matrix **classicalGramSchmidt(matrix *A)
{
    int m = A->rows;
    int n = A->cols;
    matrix **QR = malloc(sizeof(matrix *) * 2);
    QR[0] = makeMatrix(m, n);
    QR[1] = makeMatrix(n, n);
    matrix *R = QR[1];
    matrix *Q = QR[0];
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < m; i++)
        {
            Q->data[i][j] = A->data[i][j];
        }
        for (int k = 0; k < j; k++)
        {
            double dot = dotProduct(getColumn(Q,k),getColumn(Q,j));
            R->data[k][j] = dot;
            for (int i = 0; i < m; i++)
            {
                Q->data[i][j] -= dot * Q->data[i][k];
            }
        }
        double norm = dotProduct(getColumn(Q,j),getColumn(Q,j));
        norm = sqrt(norm);
        R->data[j][j] = norm;
        for (int i = 0; i < m; i++)
        {
            Q->data[i][j] /= norm;
        }
    }
    return QR;
}
/*
This function is an implementation of the modified Gram-Schmidt algorithm for QR decomposition.
It takes in a matrix A and returns a matrix Q and a matrix R such that A = QR.
It returns a pointer to an array of matrices where the first matrix is Q and the second matrix is R.
*/
matrix **modifiedGramSchmidt(matrix *A) // returning Q^T fix this
{

    int m = A->rows;
    int n = A->cols;
    matrix **QR = malloc(sizeof(matrix *) * 2);
    QR[1] = makeMatrix(n, n);
    matrix *R = QR[1];
    matrix *Q = copyMatrix(A);
    for (int i = 0; i < n; i++)
    {
        // r_ii = ||a_i||
        double r_ii = sqrt(dotProduct(getColumn(Q, i), getColumn(Q, i)));
        R->data[i][i] = r_ii;
        // q_i = a_i / r_ii
        for (int j = 0; j < m; j++)
        {
            Q->data[j][i] = Q->data[j][i] / r_ii;
        }
        for (int j = i + 1; j < n; j++)
        {
            // r_ij = q_i^T * a_j
            double r_ij = dotProduct(getColumn(Q, i), getColumn(A, j));
            R->data[i][j] = r_ij;
            // a_j = a_j - r_ij * q_i
            for (int k = 0; k < m; k++)
            {
                Q->data[k][j] = Q->data[k][j] - r_ij * Q->data[k][i];
            }
        }
    }
    QR[0] = copyMatrix(Q);
    return QR;
}
/*
This function is an implementation of the Householder algorithm for QR decomposition. It fully computes A=QR
and thus should not be used in practice
*/
matrix **inefficientHouseholder(matrix *A)
{
    // This function takes in a matrix A and returns a matrix Q and a matrix R such that A = QR
    // It returns a pointer to an array of matrices where the first matrix is Q and the second matrix is R
    matrix **QR = malloc(sizeof(matrix *) * 2);
    QR[0] = identityMatrix(A->rows, A->rows);
    QR[1] = copyMatrix(A);
    for (int k = 0; k < A->rows; k++)
    {
        vector *x = getColumn(QR[1], k);
        vector *v = makeVector(x->dimension);
        for (int i = 0; i < x->dimension; i++)
        {
            if (i < k)
            {
                v->components[i] = 0;
            }
            else
            {
                v->components[i] = x->components[i];
            }
        }
        double vNorm = sqrt(dotProduct(v, v));

        vector *e = makeVector(v->dimension);
        for (int i = 0; i < e->dimension; i++)
        {
            e->components[i] = 0;
        }
        e->components[k] = 1;

        v = vectorSubtract(v, scalarMultiplyVector(e, vNorm));
        v = makeUnitVector(v);

        matrix *H = identityMatrix(A->rows, A->rows);

        for (int i = 0; i < v->dimension; i++)
        {
            for (int j = 0; j < v->dimension; j++)
            {
                H->data[i][j] = H->data[i][j] - 2 * v->components[i] * v->components[j];
            }
        }

        QR[1] = matrixMultiply(H, QR[1]);
        QR[0] = matrixMultiply(QR[0], transpose(H));
    }
    return QR;
}

matrix **givens(matrix *A)
{
}

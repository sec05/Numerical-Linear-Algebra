#pragma once
typedef struct
{
    int dimension;
    double *components;
} vector;

vector* makeVector(int dimension);
vector* makeVectorFromData(int dimension, double* data);
vector* randomVector(int dimension);
void freeVector(vector *v);
void printVector(vector *v);
vector* scalarMultiplyVector(vector *v, double scalar);
vector* vectorAdd(vector *v1, vector *v2);
vector* vectorSubtract(vector *v1, vector *v2);
double dotProduct(vector *v1, vector *v2);
vector* makeUnitVector(vector *v);




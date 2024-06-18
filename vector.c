
#include "vector.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
vector *makeVectorFromData(int dimension, double* data){
    vector *v = malloc(sizeof(vector));
    v->dimension = dimension;
    v->components = malloc(sizeof(double) * dimension);
    for(int i = 0; i < dimension; i++){
        v->components[i] = data[i];
    }
    return v;
}
vector *makeVector(int length)
{
    vector *v = malloc(sizeof(vector));
    v->dimension = length;
    v->components = malloc(sizeof(double) * length);
    return v;
}
vector *randomVector(int length)
{
    vector *v = makeVector(length);
    for (int i = 0; i < length; i++)
    {
        v->components[i] = (double)rand() / (double)RAND_MAX;
    }
    return v;
}
void printVector(vector *v)
{
    printf("[");
    for (int i = 0; i < v->dimension; i++)
    {
        printf("%f", v->components[i]);
        if (i != v->dimension - 1)
        {
            printf(", ");
        }
    }
    printf("]\n");
}
void freeVector(vector *v)
{
    free(v->components);
    free(v);
}
vector *scalarMultiplyVector(vector *v, double scalar)
{
    vector *data = makeVector(v->dimension);
    for (int i = 0; i < v->dimension; i++)
    {
        data->components[i] = v->components[i] * scalar;
    }
    return data;
}
vector *vectorAdd(vector *v, vector *w)
{
    vector *data = makeVector(v->dimension);
    for (int i = 0; i < v->dimension; i++)
    {
        data->components[i] = v->components[i] + w->components[i];
    }
    return data;
}
vector *vectorSubtract(vector *v, vector *w)
{
    vector *data = makeVector(v->dimension);
    for (int i = 0; i < v->dimension; i++)
    {
        data->components[i] = v->components[i] - w->components[i];
    }
    return data;
}
double dotProduct(vector *v, vector *w)
{
    double sum = 0;
    int dimension = v->dimension;
    for (int i = 0; i < dimension; i++)
    {
        sum += v->components[i] * w->components[i];
    }
    return sum;
}

vector *makeUnitVector(vector *v){
    if(dotProduct(v, v) == 0){
        return v;
    }
    return scalarMultiplyVector(v, 1.0 / sqrt(dotProduct(v, v)));
}

/*
partialVectorCopy(vector *v, vector *w, int length, int index)
Parameters:
    vector *v: The vector to be copied
    vector *w: The vector to be copied to
    int length: how many elements to copy
    int index: where to start copying
*/
vector *partialVectorCopy(vector *v, vector *w, int length, int index)
{
    if(length > v->dimension - index){
        printf("Error: length is too long\n");
    }
    for (int i = 0; i < length; i++)
    {
        w->components[i] = v->components[i + index];
    }
    return w;
}

/*
partialDotProduct(vector *v, vector *w, int length, int index)
Parameters:
    vector *v: The first vector
    vector *w: The second vector
    int length: how many elements to multiply
    int index: where to start multiplying
*/
double partialDotProduct(vector *v, vector *w, int length, int index)
{
    double sum = 0;
    if(length > v->dimension - index){
        printf("Error: length is too long\n");
    }
    for (int i = 0; i < length; i++)
    {
        sum += v->components[i + index] * w->components[i + index];
    }
    return sum;
}
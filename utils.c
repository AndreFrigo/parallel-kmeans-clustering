#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>
#include <math.h>
#include <float.h>
#define EMPTY -1.0

void printMatrix(int nrow, int ncol, double *dataMatrix){
    int i,j;
    for (i = 0; i < nrow; i++){
        for (j = 0; j < ncol; j++){
            printf("%.2f\t\t", dataMatrix[i*ncol+j]);
        }
        printf("\n");
    }
    printf("----------------------------------------------------------\n");
}

void addEmptyRow(int row, int ncol, double *dataMatrix){
    int i = 0;
    for(;i<ncol+1;i++){
        dataMatrix[row*ncol+i] = EMPTY;
    }
}

bool isEmptyRow(int row, int ncol, double *dataMatrix){
    int i = 0;
    for(;i<ncol;i++){
        if(dataMatrix[row*ncol+i]!=EMPTY){
            return false;
        }
    }
    return true;
}


double distance(int r0, int r1, int ncol, double *matrix, double *centroids){
    double res = 0;
    int i;
    for(i=0;i<ncol;i++){
        res += pow((double)(matrix[r0*ncol+i] - centroids[r1*ncol+i]), 2);
    }
    return sqrt(res);
}

int chooseCluster(int p, int k, int ncol, double *matrix, double *centroids){
    int c = -1;
    double dist = DBL_MAX;
    int i;
    double d;
    //in case of non existing point return -1
    if(matrix[p*ncol+1] == EMPTY) return -1;
    
    for(i=0;i<k;i++){
        d = distance(p, i, ncol, matrix, centroids);
        if(d<dist){
            dist = d;
            c = i;
        }
    }
    return c;
}

void matrixSum(int nrow, int ncol, double *matrix1, double *matrix2){
    int i;
    for(i=0;i<nrow*ncol;i++){
        matrix1[i] += matrix2[i];
    }
}

void vectorSum(int dim, int *dest, int *source){
    int i;
    for(i=0;i<dim;i++){
        dest[i] += source[i];
    }
}


void matrixMean(int omp, int nrow, int ncol, int *vec, double *matrix){
    int i;
    #pragma omp parallel for num_threads(omp)
    for(i=0;i<nrow;i++){
        int j;
        for(j=0;j<ncol;j++){
            matrix[i*ncol+j] /= vec[i];
        } 
    }
}

bool isInArray(int elem, int dim, int *array){
    int i;
    for(i=0;i<dim;i++){
        if(array[i] == -1) return false;
        if(array[i] == elem) return true;
    }
    return false;
}

bool isSameDouble(double f0, double f1, double error){
    if(f0-f1>error || f0-f1<-1*error) return false;
    return true;
}

bool stopExecution(int omp, int nrow, int ncol, double *centroids, double *sumpoints){
    int i;
    bool ret = true;
    #pragma omp parallel for num_threads(omp)
    for(i=0;i<nrow;i++){
        int j;
        for(j=0;j<ncol;j++){
            if (!isSameDouble(centroids[i*ncol+j], sumpoints[i*ncol+j], 0.1)) ret = false;
            centroids[i*ncol+j] = sumpoints[i*ncol+j];
        }
    }
    return ret;
}
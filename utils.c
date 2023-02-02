#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>
#include <math.h>
#include <float.h>
#define EMPTY -1.0

void printMatrix(int nrow, int ncol, float *dataMatrix){
    int i,j;
    printf("PRINTING MATRIX\n");
    for (i = 0; i < nrow; i++){
        for (j = 0; j < ncol; j++){
            printf("%.2f\t\t", dataMatrix[i*ncol+j]);
        }
        printf("\n");
    }
    printf("----------------------------------------------------------\n");
}

void addEmptyRow(int row, int ncol, float *dataMatrix){
    int i = 1;
    dataMatrix[row*(ncol+1)] = -1.0;
    for(;i<ncol+1;i++){
        dataMatrix[row*(ncol+1)+i] = EMPTY;
    }
}

bool isEmptyRow(int row, int ncol, float *dataMatrix){
    int i = 0;
    for(;i<ncol;i++){
        if(dataMatrix[row*ncol+i]!=EMPTY){
            return false;
        }
    }
    return true;
}

void zeroMatrix(int omp, int nrow, int ncol, float *matrix){
    // int i,j;
    // for(i=0;i<nrow;i++){
    //     for(j=0;j<ncol;j++){
    //         matrix[i*ncol+j] = 0.0;
    //     }
    // }
    int i;
    #pragma omp parallel for num_threads(omp)
    for(i=0;i<nrow*ncol;i++){
        matrix[i] = 0.0;
    }

    return;
}

float distance(int r0, int r1, int ncol, float *matrix, float *centroids){
    double res = 0;
    int i;
    for(i=0;i<ncol;i++){
        res += pow((double)(matrix[r0*(ncol+1)+i+1] - centroids[r1*ncol+i]), 2);
    }
    return (float) sqrt(res);
}

int chooseCluster(int p, int k, int ncol, float *matrix, float *centroids){
    int c = -1;
    float dist = FLT_MAX;
    int i;
    float d;
    //in case of non existing point return -1
    if(matrix[p*(ncol+1)+1] == EMPTY) return -1;
    
    for(i=0;i<k;i++){
        d = distance(p, i, ncol, matrix, centroids);
        if(d<dist){
            dist = d;
            c = i;
        }
    }
    return c;
}

void matrixSum(int nrow, int ncol, float *matrix1, float *matrix2){
    int i;
    for(i=0;i<nrow*ncol;i++){
        matrix1[i] += matrix2[i];
    }
}

void matrixSumParallel(int omp, int nrow, int ncol, float *matrix1, float *matrix2){
    int i;
    #pragma omp parallel for num_threads(omp)
    for(i=0;i<nrow;i++){
        int j;
        for(j=0;j<ncol;j++){
            matrix1[i*ncol+j] = matrix2[i*ncol+j];
        }
    }
}

void matrixMean(int omp, int nrow, int ncol, float *matrix){
    int i;
    #pragma omp parallel for num_threads(omp)
    for(i=0;i<nrow;i++){
        int j;
        for(j=1;j<ncol;j++) matrix[i*ncol+j] /= matrix[i*ncol];
        matrix[i*ncol] = 1.0; 
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

bool isSameFloat(float f0, float f1, float error){
    if(f0-f1>error || f0-f1<-1*error) return false;
    return true;
}

bool stopExecution(int omp, int nrow, int ncol, float *centroids, float *sumpoints){
    int i;
    bool ret = true;
    #pragma omp parallel for num_threads(omp)
    for(i=0;i<nrow;i++){
        int j;
        for(j=0;j<ncol;j++){
            if (!isSameFloat(centroids[i*ncol+j], sumpoints[i*(ncol+1)+j+1], 0.1)) ret = false;
            centroids[i*ncol+j] = sumpoints[i*(ncol+1)+j+1];
        }
    }
    return ret;
}
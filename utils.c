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
    int i = 0;
    for(;i<ncol;i++){
        dataMatrix[row*ncol+i] = EMPTY;
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

float distance(int r0, int r1, int ncol, float *matrix0, float *matrix1){
    double res = 0;
    int i;
    for(i=0;i<ncol;i++){
        res += pow((double)(matrix0[r0*ncol+i] - matrix1[r1*ncol+i]), 2);
    }
    return (float) sqrt(res);
}

int chooseCluster(int p, int k, int ncol, float *matrix, float *centroids){
    int c = -1;
    float dist = FLT_MAX;
    int i;
    float d;
    //in case of non existing point return -1
    if(matrix[p*ncol] == EMPTY) return -1;
    
    for(i=0;i<k;i++){
        d = distance(p, i, ncol, matrix, centroids);
        if(d<dist){
            dist = d;
            c = i;
        }
    }
    return c;
}
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
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
#include <stdio.h>
#include <stdlib.h>

void printMatrix(int nrow, int ncol, float* dataMatrix){
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
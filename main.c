#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "functions.h"

int main(int argc, char *argv[]){
    //process id
    int my_rank;
    //number of processes
    int n_proc;
    
    int nrow, ncol;

    //2 arguments are expected: first the filename of the dataset, second the number of MPI processes
    if (argc < 3){
        printf("Error! Wrong number of arguments");
        return -1;
    }
    int nproc = atoi(argv[2]);

    char *filename = argv[1];

    nrow = getRows(filename);
    ncol = getCols(filename);
    printf("Num rows: %d\nNum cols: %d\n", nrow, ncol);

    //Dataset matrix, usage:
    //dataMatrix[i * ncol + j] corresponds to dataMatrix[i][j]
    //dataMatrix[i + nrow * j] corresponds to dataMatrix[i][j]
    float* dataMatrix=NULL;
    dataMatrix = (float *)malloc(nrow * ncol * sizeof(float));

    readFile(filename, nrow, ncol, dataMatrix);

    printMatrix(nrow, ncol, dataMatrix);


}

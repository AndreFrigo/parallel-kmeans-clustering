#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "functions.h"

int main(int argc, char *argv[]){
    //process id
    int my_rank;
    //number of processes
    int n_proc;


    //Dataset matrix, usage:
    //dataMatrix[i * ncol + j] corresponds to dataMatrix[i][j]
    //dataMatrix[i + nrow * j] corresponds to dataMatrix[i][j]

    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD , &n_proc);

    int nrow, ncol;
    float *dataMatrix=NULL;

    if(my_rank == 0){
        //2 arguments are expected: first the filename of the dataset, second the number of MPI processes
        if (argc < 3){
            printf("Error! Wrong number of arguments");
            return -1;
        }
        int nproc = atoi(argv[2]);

        char *filename = argv[1];

        printf("P0\n");
        nrow = getRows(filename);
        ncol = getCols(filename);
        //make the number of rows divisible by the number of MPI processes used
        if(nrow%n_proc!=0){
            nrow += n_proc-nrow%n_proc;
        }

        //broadcast nrow and ncol
        MPI_Bcast(&nrow, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&ncol, 1, MPI_INT, 0, MPI_COMM_WORLD);

        dataMatrix = (float *)malloc(nrow * ncol * sizeof(float));
        
        readFile(filename, nrow, ncol, dataMatrix);

        // printMatrix(nrow, ncol, dataMatrix);

    }

    //broadcast nrow and ncol
    MPI_Bcast(&nrow, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ncol, 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("PROC %d\nNum rows: %d\nNum cols: %d\n", my_rank, nrow, ncol);

    //scatter the dataMatrix
    int scatterRow=nrow/n_proc;
    float* recvMatrix = (float *)malloc(scatterRow * ncol * sizeof(float));
    //scatter the matrix
    MPI_Scatter(dataMatrix, scatterRow*ncol, MPI_FLOAT, recvMatrix, scatterRow*ncol, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // printf("PROCESS %d\n", my_rank);
    // printMatrix(scatterRow, ncol, recvMatrix);

    if(my_rank==0) free(dataMatrix);


    return 0;
}

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
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

    int nrowold, nrow, ncol, k;
    float *dataMatrix=NULL;
    float *centroids=NULL;

    if(my_rank == 0){
        //2 arguments are expected: first the filename of the dataset, second the number of MPI processes
        if (argc < 4){
            printf("Error! Wrong number of arguments");
            return -1;
        }
        char *filename = argv[1];
        int nproc = atoi(argv[2]);
        k = atoi(argv[3]);
        


        printf("P0\n");
        printf("K: %d\n", k);
        nrowold = getRows(filename);
        ncol = getCols(filename);
        //make the number of rows divisible by the number of MPI processes used
        if(nrowold%n_proc!=0){
            nrow = nrowold + n_proc-nrowold%n_proc;
        }

        // //broadcast nrow and ncol
        // MPI_Bcast(&nrow, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // MPI_Bcast(&ncol, 1, MPI_INT, 0, MPI_COMM_WORLD);

        dataMatrix = (float *)malloc(nrow * ncol * sizeof(float));
        
        
        readFile(filename, nrow, ncol, dataMatrix);

        // printMatrix(nrow, ncol, dataMatrix);

    }

    //broadcast nrow and ncol
    MPI_Bcast(&nrow, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ncol, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("PROC %d\nNum rows: %d\nNum cols: %d\nK: %d\n", my_rank, nrow, ncol, k);

    //scatter the dataMatrix
    int scatterRow=nrow/n_proc;
    float* recvMatrix = (float *)malloc(scatterRow * ncol * sizeof(float));
    //scatter the matrix
    MPI_Scatter(dataMatrix, scatterRow*ncol, MPI_FLOAT, recvMatrix, scatterRow*ncol, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // printf("PROCESS %d\n", my_rank);
    // printMatrix(scatterRow, ncol, recvMatrix);

    centroids = (float *)malloc(k * ncol * sizeof(float));
    if(my_rank==0){
        printf("INSIDE CRITICAL SECTION, K: %d\n", k);
        //choose randomly k centroids
        int i;
        srand(time(NULL));
        for(i=0;i<k;i++){
            int r = rand()%nrowold;
            int c;
            printf("Random: row: %d\n", r);
            for(c=0;c<ncol;c++){
                centroids[i*ncol+c] = dataMatrix[r*ncol+c];
            }
        }
        // printMatrix(k, ncol, centroids);
        free(dataMatrix);
    }
    MPI_Bcast(centroids, k*ncol, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // printf("PROC %d\n", my_rank);
    // printMatrix(k, ncol, centroids);


    return 0;
}

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <stdbool.h>
#include "functions.h"

int main(int argc, char *argv[]){
    //process id
    int my_rank;
    //number of processes
    int n_proc;
    //number of threads
    int omp;

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
        omp = atoi(argv[2]);
        k = atoi(argv[3]);
        


        // printf("P0\n");
        // printf("K: %d\n", k);
        //useful to know how many points I have
        nrowold = getRows(filename);
        // printf("Num real rows: %d\n", nrowold);
        ncol = getCols(filename);
        //make the number of rows divisible by the number of MPI processes used
        if(nrowold%n_proc!=0){
            nrow = nrowold + n_proc-nrowold%n_proc;
        }else{
            nrow = nrowold;
        }

        // //broadcast nrow and ncol
        // MPI_Bcast(&nrow, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // MPI_Bcast(&ncol, 1, MPI_INT, 0, MPI_COMM_WORLD);

        dataMatrix = (float *)malloc(nrow * (ncol+1) * sizeof(float));
        
        
        readFile(filename, nrow, ncol, dataMatrix);

        // printMatrix(nrow, ncol+1, dataMatrix, false);

    }

    //broadcast nrow and ncol
    MPI_Bcast(&nrow, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ncol, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&omp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // printf("PROC %d\nNum rows: %d\nNum cols: %d\nK: %d\n", my_rank, nrow, ncol, k);

    //scatter the dataMatrix
    int scatterRow=nrow/n_proc;
    float* recvMatrix = (float *)malloc(scatterRow * (ncol+1) * sizeof(float));
    //scatter the matrix
    MPI_Scatter(dataMatrix, scatterRow*(ncol+1), MPI_FLOAT, recvMatrix, scatterRow*(ncol+1), MPI_FLOAT, 0, MPI_COMM_WORLD);

    // printf("PROCESS %d\n", my_rank);
    // printMatrix(scatterRow, ncol+1, recvMatrix, false);

    centroids = (float *)malloc(k * ncol * sizeof(float));
    //for each process a matrix that stores the sum of all points of each cluster and the number of points
    float *sumpoints = NULL;
    //same as before, but for P0 to store the sum of all the sumpoints matrices and compare the results with the centroids of the previous iteration
    float *sumpointsP0 = NULL;
    if(my_rank==0){
        //choose randomly k centroids
        int i;
        //TODO
        // srand(time(NULL));
        srand(1);
        int alreadySelected[k];
        for (i=0;i<k;i++) alreadySelected[i] = -1;
        for(i=0;i<k;i++){
            int r = rand()%nrowold;
            if(!isInArray(r, k, alreadySelected)){
                int c;
                // printf("Random: row: %d\n", r);
                for(c=0;c<ncol;c++){
                    centroids[i*ncol+c] = dataMatrix[r*(ncol+1)+c+1];
                }
                alreadySelected[i] = r;
            }else{
                i--;
            }
        }
        // printMatrix(k, ncol, centroids, true);
        free(dataMatrix);
        sumpointsP0 = (float *)malloc(k * (ncol+1) * sizeof(float));
    }

    // //variable to check when to stop the algorithm
    // bool stop = false;

    // while(!stop){

    // }

    MPI_Bcast(centroids, k*ncol, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // printf("PROC %d\n", my_rank);
    // printMatrix(k, ncol, centroids, true);

    
    sumpoints = (float *)malloc(k * (ncol+1) * sizeof(float));

    zeroMatrix(omp, k, ncol+1, sumpoints);
    //TODO: handle clusters structure (keep the points divided by clusters to return the final grouping???)
    #pragma omp parallel num_threads(omp)
    {
        int i;
        //matrix to store sums for each threads (scope private)
        float partialMatrix[k*(ncol+1)];
        for(i=0;i<k*(ncol+1);i++) partialMatrix[i] = 0;
        int res;
        //store partial sums in the private matrix
        #pragma omp for schedule(static, (int) scatterRow/omp)
        for(i=0; i<scatterRow;i++){
            res = chooseCluster(i, k, ncol, recvMatrix, centroids);
            if(res!=-1){
                int j;
                partialMatrix[res*(ncol+1)]++;
                for(j=1;j<ncol+1;j++){
                    partialMatrix[res*(ncol+1)+j] += recvMatrix[i*(ncol+1)+j];
                    recvMatrix[i*(ncol+1)] = res;
                }
            }
            // printf("Proc %d Thread %d Point %d: cluster %d\n", my_rank, omp_get_thread_num(), i, chooseCluster(i, k, ncol, recvMatrix, centroids));
        }
        //sum up the different private matrices in the shared one
        #pragma omp critical
        matrixSum(k, ncol+1, sumpoints, partialMatrix);
    }
    // printf("PROC %d\n",my_rank);
    // printMatrix(k, ncol+1, sumpoints, true);

    MPI_Reduce(sumpoints, sumpointsP0, k*(ncol+1), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(my_rank==0){
        // printMatrix(k, ncol+1, sumpointsP0, true);
        matrixMean(omp, k, ncol+1, sumpointsP0);
        printf("Previous centroids\n");
        printMatrix(k, ncol, centroids, true);
        int i,j;
        for(i=0; i<k; i++){
            for(j=0;j<ncol;j++) centroids[i*ncol+j] = sumpointsP0[i*(ncol+1)+j+1];
        }
        printf("Actual centroids\n");
        printMatrix(k, ncol, centroids, true);
    } 
    return 0;
}

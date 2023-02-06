#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <stdbool.h>
#include <sys/time.h>
#include "functions.h"
#define DEBUG 1
#define NUMITER 0

int main(int argc, char *argv[]){
    //process id
    int my_rank;
    //number of processes
    int n_proc;
    //number of threads
    int omp;

    //Dataset matrix, usage:
    //dataMatrix[i * ncol + j] corresponds to dataMatrix[i][j]

    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD , &n_proc);

    int nrowold, nrow, ncol, k, cont;
    float *dataMatrix=NULL;
    float *centroids=NULL;
    struct timeval start, afterReading, beforeWhile, end;
    int *mapping=NULL;
    int *counter=NULL;
    int *counterP0 = NULL;

    if(my_rank == 0){
        //3 arguments are expected: first the filename of the dataset, second the number of OMP processes, third the number of clusters
        if (argc < 4){
            printf("Error! Wrong number of arguments");
            return -1;
        }
        char *filename = argv[1];
        omp = atoi(argv[2]);
        k = atoi(argv[3]);
        cont = 0;

        gettimeofday(&start, NULL);

        nrowold = getRows(filename);
        ncol = getCols(filename);
        //make the number of rows divisible by the number of MPI processes used
        if(nrowold%n_proc!=0){
            nrow = nrowold + n_proc-nrowold%n_proc;
        }else{
            nrow = nrowold;
        }

        dataMatrix = (float *)malloc(nrow * ncol * sizeof(float));
        readFile(filename, nrow, ncol, dataMatrix);
        mapping = (int *)malloc(nrow * sizeof(int));
        // int i;
        // for(i=0;i<nrow;i++) mapping[i] = -1;
        //save start time after reading the dataset
        gettimeofday(&afterReading, NULL);
    }

    //broadcast nrow and ncol
    MPI_Bcast(&nrow, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ncol, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&omp, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //scatter the dataMatrix
    int scatterRow=nrow/n_proc;
    //scatter the matrix and the mapping vector
    float* recvMatrix = (float *)malloc(scatterRow * ncol * sizeof(float));
    MPI_Scatter(dataMatrix, scatterRow*ncol, MPI_FLOAT, recvMatrix, scatterRow*ncol, MPI_FLOAT, 0, MPI_COMM_WORLD);
    int *recvMapping = (int *)malloc(scatterRow * sizeof(int));
    MPI_Scatter(mapping, scatterRow, MPI_INT, recvMapping, scatterRow, MPI_INT, 0, MPI_COMM_WORLD);

    centroids = (float *)malloc(k * ncol * sizeof(float));
    //for each process a matrix that stores the sum of all points of each cluster and the number of points
    float *sumpoints = NULL;
    //same as before, but for P0 to store the sum of all the sumpoints matrices and compare the results with the centroids of the previous iteration
    float *sumpointsP0 = NULL;
    if(my_rank==0){
        //choose randomly k centroids
        int i;
        if(DEBUG){
            srand(1);
        }else{
            srand(time(NULL));
        }
        int alreadySelected[k];
        for (i=0;i<k;i++) alreadySelected[i] = -1;
        for(i=0;i<k;i++){
            int r = rand()%nrowold;
            if(!isInArray(r, k, alreadySelected)){
                int c;
                for(c=0;c<ncol;c++){
                    centroids[i*ncol+c] = dataMatrix[r*ncol+c];
                }
                alreadySelected[i] = r;
            }else{
                i--;
            }
        }
        if(DEBUG){
            printf("Generated centroids:\n");
            printMatrix(k, ncol, centroids);
        }
        sumpointsP0 = (float *)malloc(k * ncol * sizeof(float));
        counterP0 = (int *)malloc(k * sizeof(int));
    }

    sumpoints = (float *)malloc(k * ncol * sizeof(float));
    counter = (int *)malloc(k * sizeof(int));
    if(my_rank==0) gettimeofday(&beforeWhile, NULL);
    //variable to check when to stop the algorithm
    bool stop = false;
    while(!stop){
        MPI_Bcast(centroids, k*ncol, MPI_FLOAT, 0, MPI_COMM_WORLD);
        #pragma omp parallel num_threads(omp) shared(sumpoints)
        {
            int i;
            //reset the sumpoints matrix
            #pragma omp for
            for(i=0;i<k*ncol;i++) sumpoints[i] = 0.0;
            #pragma omp for
            for(i=0;i<k;i++) counter[i] = 0;
            //matrix to store sums for each threads (scope private)
            float partialMatrix[k*ncol];
            int partialCounter[k];
            for(i=0;i<k*ncol;i++) partialMatrix[i] = 0.0;
            for(i=0;i<k;i++) partialCounter[i] = 0;
            int res;
            //store partial sums in the private matrix
            #pragma omp for nowait schedule(static, 1)
            for(i=0; i<scatterRow;i++){
                res = chooseCluster(i, k, ncol, recvMatrix, centroids);
                if(res!=-1){
                    int j;
                    partialCounter[res]++;
                    recvMapping[i] = res;
                    //TODO: test with parallel for
                    for(j=0;j<ncol;j++){
                        partialMatrix[res*ncol+j] += recvMatrix[i*ncol+j];
                    }
                }
            }
            //sum up the different private matrices in the shared one
            #pragma omp critical
            {
                matrixSum(k, ncol, sumpoints, partialMatrix);
                vectorSum(k, counter, partialCounter);
            }

        }

        MPI_Reduce(sumpoints, sumpointsP0, k*ncol, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(counter, counterP0, k, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if(my_rank==0){
            cont++;
            matrixMean(omp, k, ncol, counterP0, sumpointsP0);
            if (stopExecution(omp, k, ncol, centroids, sumpointsP0) || (NUMITER > 0 && cont >= NUMITER)){
                stop = true;
                printf("Stop execution, printing final centroids\n");
                printMatrix(k, ncol, centroids);
            }else{
                // for debugging iterations
                printf("Iteration %d, actual centroids:\n", cont);
                printMatrix(k, ncol, centroids);
            }
        } 
        //broadcast the stopping condition, all processes have to stop in the same cycle
        MPI_Bcast(&stop, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

    }
    free(sumpoints);
    free(sumpointsP0);
    free(counter);
    free(counterP0);
    
    MPI_Gather(recvMatrix, scatterRow*ncol, MPI_FLOAT, dataMatrix, scatterRow*ncol, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(recvMapping, scatterRow, MPI_FLOAT, mapping, scatterRow, MPI_FLOAT, 0, MPI_COMM_WORLD);
    free(recvMatrix);
    free(recvMapping);
    if(my_rank==0){
        //calculate execution time before printing the matrix, but after gathering it
        gettimeofday(&end, NULL);

        //print statistics
        printf("NROW: %d, NCOL: %d, NITER: %d\n", nrowold, ncol, cont);
        printf("EXECUTION TIME FOR DIFFERENT PHASES IN MICROSECONDS\n");
        printf("TOTAL EXECUTION TIME: %ld\n", ((end.tv_sec*1000000 + end.tv_usec) -(start.tv_sec*1000000 + start.tv_usec)));
        printf("READING DATASET TIME: %ld\n", ((afterReading.tv_sec*1000000 + afterReading.tv_usec) -(start.tv_sec*1000000 + start.tv_usec)));
        printf("AVERAGE CYCLIC EXECUTION TIME: %ld\n", ((end.tv_sec*1000000 + end.tv_usec) -(beforeWhile.tv_sec*1000000 + beforeWhile.tv_usec))/cont);
        // no need to remove fake points from the matrix, the gather works by rank order, so it's enough to use nrowold instead of nrow
        // printMatrix(nrowold, ncol, dataMatrix);
    }
    MPI_Finalize();
    return 0;
}

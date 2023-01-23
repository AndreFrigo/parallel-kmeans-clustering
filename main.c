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

    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD , &n_proc);

    int nrowold, nrow, ncol, k;
    float *dataMatrix=NULL;
    float *centroids=NULL;

    if(my_rank == 0){
        //3 arguments are expected: first the filename of the dataset, second the number of OMP processes, third the number of clusters
        if (argc < 4){
            printf("Error! Wrong number of arguments");
            return -1;
        }
        char *filename = argv[1];
        omp = atoi(argv[2]);
        k = atoi(argv[3]);
        
        nrowold = getRows(filename);
        ncol = getCols(filename);
        //make the number of rows divisible by the number of MPI processes used
        if(nrowold%n_proc!=0){
            nrow = nrowold + n_proc-nrowold%n_proc;
        }else{
            nrow = nrowold;
        }

        dataMatrix = (float *)malloc(nrow * (ncol+1) * sizeof(float));
        readFile(filename, nrow, ncol, dataMatrix);

    }

    //broadcast nrow and ncol
    MPI_Bcast(&nrow, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ncol, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&omp, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //scatter the dataMatrix
    int scatterRow=nrow/n_proc;
    float* recvMatrix = (float *)malloc(scatterRow * (ncol+1) * sizeof(float));
    //scatter the matrix
    MPI_Scatter(dataMatrix, scatterRow*(ncol+1), MPI_FLOAT, recvMatrix, scatterRow*(ncol+1), MPI_FLOAT, 0, MPI_COMM_WORLD);


    centroids = (float *)malloc(k * ncol * sizeof(float));
    //for each process a matrix that stores the sum of all points of each cluster and the number of points
    float *sumpoints = NULL;
    //same as before, but for P0 to store the sum of all the sumpoints matrices and compare the results with the centroids of the previous iteration
    float *sumpointsP0 = NULL;
    if(my_rank==0){
        //choose randomly k centroids
        int i;
        srand(time(NULL));
        // srand(1);
        int alreadySelected[k];
        for (i=0;i<k;i++) alreadySelected[i] = -1;
        for(i=0;i<k;i++){
            int r = rand()%nrowold;
            if(!isInArray(r, k, alreadySelected)){
                int c;
                for(c=0;c<ncol;c++){
                    centroids[i*ncol+c] = dataMatrix[r*(ncol+1)+c+1];
                }
                alreadySelected[i] = r;
            }else{
                i--;
            }
        }
        free(dataMatrix);
        sumpointsP0 = (float *)malloc(k * (ncol+1) * sizeof(float));
    }

    sumpoints = (float *)malloc(k * (ncol+1) * sizeof(float));
    //variable to check when to stop the algorithm
    bool stop = false;
    while(!stop){
        MPI_Bcast(centroids, k*ncol, MPI_FLOAT, 0, MPI_COMM_WORLD);

        zeroMatrix(omp, k, ncol+1, sumpoints);

        #pragma omp parallel num_threads(omp)
        {
            int i;
            //matrix to store sums for each threads (scope private)
            float partialMatrix[k*(ncol+1)];
            for(i=0;i<k*(ncol+1);i++) partialMatrix[i] = 0.0;
            int res;
            //store partial sums in the private matrix
            #pragma omp for schedule(static, (int) scatterRow/omp)
            for(i=0; i<scatterRow;i++){
                res = chooseCluster(i, k, ncol, recvMatrix, centroids);
                if(res!=-1){
                    int j;
                    partialMatrix[res*(ncol+1)]++;
                    recvMatrix[i*(ncol+1)] = res;
                    for(j=1;j<ncol+1;j++){
                        partialMatrix[res*(ncol+1)+j] += recvMatrix[i*(ncol+1)+j];
                    }
                }
            }
            //sum up the different private matrices in the shared one
            #pragma omp critical
            matrixSum(k, ncol+1, sumpoints, partialMatrix);

        }

        MPI_Reduce(sumpoints, sumpointsP0, k*(ncol+1), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        if(my_rank==0){
            matrixMean(omp, k, ncol+1, sumpointsP0);
            if (stopExecution(omp, k, ncol, centroids, sumpointsP0)){
                stop = true;
                printf("Stop execution, printing final centroids\n");
                printMatrix(k, ncol, centroids, true);
            }else{
                // for debugging iterations
                // printf("Actual centroids\n");
                // printMatrix(k, ncol, centroids, true);
            }
        } 
        //broadcast the stopping condition, all processes have to stop in the same cycle
        MPI_Bcast(&stop, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

    }
    //out of loop, the algorithm has already finished, send the data with the choice of the cluster back to P0
    if(my_rank==0){
        dataMatrix = (float *)malloc(nrow * (ncol+1) * sizeof(float));
    }else{
        dataMatrix = NULL;
    }
    MPI_Gather(recvMatrix, scatterRow*(ncol+1), MPI_FLOAT, dataMatrix, scatterRow*(ncol+1), MPI_FLOAT, 0, MPI_COMM_WORLD);
    free(recvMatrix);
    if(my_rank==0){
        printMatrix(nrow, ncol+1, dataMatrix, true);
    }
    MPI_Finalize();
    return 0;
}

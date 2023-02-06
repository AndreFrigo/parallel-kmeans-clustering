#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <sys/time.h>
#include <math.h>
#include <float.h>
#include <string.h>
#define EMPTY -1.0
#define MAXCHAR 10000
#define DEBUG 1
#define NUMITER 0


int getRows(char *filename){
    FILE *fp;
    int count = 0;
    char c;

    fp = fopen(filename, "r");
    if (fp == NULL){
        printf("Error while opening file %s", filename);
        return 0;
    }

    for (c = getc(fp); c != EOF; c = getc(fp))
        if (c == '\n') count += 1;
 
    fclose(fp);
    //first row is the header
    return count-1;
}

int getCols(char *filename){
    FILE *fp;
    int count = 0;
    char c;

    fp = fopen(filename, "r");
    if (fp == NULL){
        printf("Error while opening file %s", filename);
        return 0;
    }

    for (c = getc(fp); c != '\n'; c = getc(fp))
        if (c == ',') count += 1;
 
    fclose(fp);
    //last element has no ',' but a '\n'
    return count+1;
}

void readFile(char *filename, int nrow, int ncol, double *dataMatrix){
    FILE *fp;
    char row[MAXCHAR];
    char *token;

    fp = fopen(filename, "r");
    if (fp == NULL){
        printf("Error while opening file %s", filename);
        return;
    }

    int r=0;
    //header row
    fgets(row, MAXCHAR, fp);

    while (feof(fp) != true){
        fgets(row, MAXCHAR, fp);
        token = strtok(row, ",");
        int c = 0;
        while(token != NULL){
            dataMatrix[r*ncol+c] = atof(token);
            token = strtok(NULL, ",");
            c++;
        }
        r++;
    }
    fclose(fp);
    //the last line of input file is an empty one
    r--;

}

double distance(int r0, int r1, int ncol, double *matrix, double *centroids){
    double res = 0;
    int i;
    for(i=0;i<ncol;i++){
        res += pow((double)(matrix[r0*ncol+i] - centroids[r1*ncol+i]), 2);
    }
    return sqrt(res);
}

int chooseCluster(int p, int k, int ncol, double *matrix, double *centroids){
    int c = -1;
    double dist = DBL_MAX;
    int i;
    double d;
    //in case of non existing point return -1
    if(matrix[p*ncol+1] == EMPTY) return -1;
    
    for(i=0;i<k;i++){
        d = distance(p, i, ncol, matrix, centroids);
        if(d<dist){
            dist = d;
            c = i;
        }
    }
    return c;
}

void zeroMatrix(int nrow, int ncol, double *matrix){
    int i;
    for(i=0;i<nrow*ncol;i++){
        matrix[i] = 0.0;
    }

    return;
}

void matrixMean(int nrow, int ncol, int *vec, double *matrix){
    int i;
    for(i=0;i<nrow;i++){
        int j;
        for(j=0;j<ncol;j++) matrix[i*ncol+j] /= vec[i];
    }
}

bool isSameDouble(double f0, double f1, double error){
    if(f0-f1>error || f0-f1<-1*error) return false;
    return true;
}

bool isInArray(int elem, int dim, int *array){
    int i;
    for(i=0;i<dim;i++){
        if(array[i] == -1) return false;
        if(array[i] == elem) return true;
    }
    return false;
}

bool stopExecution(int nrow, int ncol, double *centroids, double *sumpoints){
    int i;
    bool ret = true;
    for(i=0;i<nrow;i++){
        int j;
        for(j=0;j<ncol;j++){
            if (!isSameDouble(centroids[i*ncol+j], sumpoints[i*ncol+j], 0.1)) ret = false;
            centroids[i*ncol+j] = sumpoints[i*ncol+j];
        }
    }
    return ret;
}

void printMatrix(int nrow, int ncol, double *dataMatrix){
    int i,j;
    for (i = 0; i < nrow; i++){
        j=0;
        for (; j < ncol; j++){
            printf("%.2f\t\t", dataMatrix[i*ncol+j]);
        }
        printf("\n");
    }
    printf("----------------------------------------------------------\n");
}

int main(int argc, char *argv[]){

    int nrow, ncol, k;
    double *dataMatrix=NULL;
    double *sumpoints=NULL;
    double *centroids=NULL;
    struct timeval start, afterReading, afterRandomCentroids, beforeWhile, end;
    int *mapping=NULL;
    int *counter=NULL;

    if (argc < 3){
        printf("Error! Wrong number of arguments");
        return -1;
    }

    char *filename = argv[1];
    k = atoi(argv[2]);
    //save start time after reading the dataset
    gettimeofday(&start, NULL);

    nrow = getRows(filename);
    ncol = getCols(filename);
    dataMatrix = (double *)malloc(nrow * ncol * sizeof(double));
    mapping = (int *)malloc(nrow * sizeof(int));
    readFile(filename, nrow, ncol, dataMatrix);
    //save start time after reading the dataset
    gettimeofday(&afterReading, NULL);

    centroids = (double *)malloc(k * ncol * sizeof(double));
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
        int r = rand()%nrow;
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
    gettimeofday(&afterRandomCentroids, NULL);
    if(DEBUG){
        printf("Generated centroids:\n");
        printMatrix(k, ncol, centroids);
    }
    sumpoints = (double *)malloc(k * ncol * sizeof(double));
    counter = (int *)malloc(k * sizeof(int));
    //variable to check when to stop the algorithm
    bool stop = false;
    int cont = 0;
    int res;
    gettimeofday(&beforeWhile, NULL);
    while(!stop){
        zeroMatrix(k, ncol, sumpoints);
        int i;
        for(i=0;i<k;i++) counter[i] = 0;
        for(i=0;i<nrow;i++){
            res = chooseCluster(i, k, ncol, dataMatrix, centroids);
            if(res!=-1){
                int j;
                counter[res]++;
                mapping[i] = res;
                for(j=0;j<ncol;j++){
                    sumpoints[res*ncol+j] += dataMatrix[i*ncol+j];
                }
            }
        }

        matrixMean(k, ncol, counter, sumpoints);
        cont++;

    

        
        if (stopExecution(k, ncol, centroids, sumpoints) || (NUMITER>0 && cont>=NUMITER)){
            stop = true;
            printf("Stop execution after %d cycles, printing final centroids\n", cont);
            printMatrix(k, ncol, centroids);
        }
        // else{
        //     if(DEBUG){
        //         printf("Iteration %d, actual centroids:\n", cont);
        //         printMatrix(k, ncol, centroids);
        //     }
        // }
        
    }

    gettimeofday(&end, NULL);

    //print statistics
    printf("NROW: %d, NCOL: %d, NITER: %d\n", nrow, ncol, cont);
    printf("EXECUTION TIME FOR DIFFERENT PHASES IN MICROSECONDS\n");
    printf("TOTAL EXECUTION TIME: %ld\n", ((end.tv_sec*1000000 + end.tv_usec) -(start.tv_sec*1000000 + start.tv_usec)));
    printf("READING DATASET TIME: %ld\n", ((afterReading.tv_sec*1000000 + afterReading.tv_usec) -(start.tv_sec*1000000 + start.tv_usec)));
    printf("GENERATING RANDOM CENTROIDS TIME: %ld\n", ((afterRandomCentroids.tv_sec*1000000 + afterRandomCentroids.tv_usec) -(afterReading.tv_sec*1000000 + afterReading.tv_usec)));
    printf("AVERAGE CYCLIC EXECUTION TIME: %ld\n", ((end.tv_sec*1000000 + end.tv_usec) -(beforeWhile.tv_sec*1000000 + beforeWhile.tv_usec))/cont);

    return 0;

}
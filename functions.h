#include <stdbool.h>

//inputfunctions.c

//get number of entries
int getRows(char *filename);
//get number of attributes of each entry
int getCols(char *filename);
//read the dataset and store it in a matrix
void readFile(char *filename, int nrow, int ncol, float *dataMatrix);

//utils.c 

//print a float matrix
void printMatrix(int nrow, int ncol, float *dataMatrix);
//add an empty row to the matrix
void addEmptyRow(int row, int ncol, float *dataMatrix);
//check if a row is not significant, return true if it is empty, otherwise false
bool isEmptyRow(int row, int ncol, float *dataMatrix);
//initialize a matrix inserting 0.0 in each cell
void zeroMatrix(int omp, int nrow, int ncol, float *matrix);
//calculate distance between two points (one point of the matrix and one centroid)
float distance(int r0, int r1, int ncol, float *matrix, float *centroids);
//returns the index of the centroid nearest to a specific point, returns -1 in case of not real point
int chooseCluster(int p, int k, int ncol, float *matrix, float *centroids);
//sum two matrixes (serial)
void matrixSum(int nrow, int ncol, float *matrix1, float *matrix2);
//sum two matrixes (parallel with OMP)
void matrixSumParallel(int omp, int nrow, int ncol, float *matrix1, float *matrix2);
//calculates the weighted mean of a matrix that has as first column the weight
void matrixMean(int omp, int nrow, int ncol, float *matrix);
//check if an integer element is in an integer array
bool isInArray(int elem, int dim, int *array);
//check if two floating point values are the same
bool isSameFloat(float f0, float f1, float error);
//compare the old and new centroids to decide whether to stop
bool stopExecution(int nrow, int ncol, float *matrix0, float *matrix1);
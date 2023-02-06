#include <stdbool.h>

//inputfunctions.c

//get number of entries
int getRows(char *filename);
//get number of attributes of each entry
int getCols(char *filename);
//read the dataset and store it in a matrix
void readFile(char *filename, int nrow, int ncol, double *dataMatrix);

//utils.c 

//print a double matrix
void printMatrix(int nrow, int ncol, double *dataMatrix);
//add an empty row to the matrix
void addEmptyRow(int row, int ncol, double *dataMatrix);
//check if a row is not significant, return true if it is empty, otherwise false
bool isEmptyRow(int row, int ncol, double *dataMatrix);
//calculate distance between two points (one point of the matrix and one centroid)
double distance(int r0, int r1, int ncol, double *matrix, double *centroids);
//returns the index of the centroid nearest to a specific point, returns -1 in case of not real point
int chooseCluster(int p, int k, int ncol, double *matrix, double *centroids);
//sum two matrixes (serial)
void matrixSum(int nrow, int ncol, double *matrix1, double *matrix2);
//sum two vectors (serial)
void vectorSum(int dim, int *dest, int *source);
//calculates the weighted mean of a matrix that has as first column the weight
void matrixMean(int omp, int nrow, int ncol, int *vec, double *matrix);
//check if an integer element is in an integer array
bool isInArray(int elem, int dim, int *array);
//check if two double point values are the same
bool isSameDouble(double f0, double f1, double error);
//compare the old and new centroids to decide whether to stop (the new centroids in form of sum matrix), implemented using OMP only in P0, it also updates the current centroids
bool stopExecution(int omp, int nrow, int ncol, double *centroids, double *summatrix);
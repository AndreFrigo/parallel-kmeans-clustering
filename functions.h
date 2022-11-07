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
//calculate distance between two points (rows) of a matrix using Euclidean distance
float distance(int r0, int r1, int ncol, float *matrix0, float *matrix1);
//returns the index of the centroid nearest to a specific point
int chooseCluster(int p, int k, int ncol, float *matrix, float *centroids);
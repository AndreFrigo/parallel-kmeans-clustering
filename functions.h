//inputfunctions.c

//get number of entries
int getRows(char *filename);
//get number of attributes of each entry
int getCols(char *filename);
//read the dataset and store it in a matrix
void readFile(char * filename, int nrow, int ncol, float* dataMatrix);

//utils.c 

//print a float matrix
void printMatrix(int nrow, int ncol, float* dataMatrix);
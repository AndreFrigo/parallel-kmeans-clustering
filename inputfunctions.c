#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include "functions.h"
#define MAXCHAR 10000

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

void readFile(char *filename, int nrow, int ncol, float *dataMatrix){
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
    //add empty rows to complete the matrix
    for(;r<nrow;r++){
        addEmptyRow(r, ncol, dataMatrix);
    }

}
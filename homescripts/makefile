SHELL:=/bin/bash

parallel: parallel-kmeans-clustering/inputfunctions.c parallel-kmeans-clustering/utils.c parallel-kmeans-clustering/main.c
	module load mpich-3.2; \
	mpicc -g -Wall -fopenmp -o executable parallel-kmeans-clustering/inputfunctions.c parallel-kmeans-clustering/utils.c parallel-kmeans-clustering/main.c -lm
serial: parallel-kmeans-clustering/serial.c
	gcc -o executable parallel-kmeans-clustering/serial.c -lm
clean:
	rm -f script.s*
	rm -f executable

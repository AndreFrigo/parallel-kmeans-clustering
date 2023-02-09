IN ORDER TO RUN THE SCRIPT, THE FILES ON THIS FOLDER HAS TO BE PUT IN HOME

To run the program:
1) 'make' or 'make parallel' to compile and create executable ('make serial' for the serial implementation)
2) ./parallel.sh [nodes] [cpus/node] [dataset] [k] (./serial.sh [dataset] [k] for serial implementation)

To read the results: 'cat script.sh.oXXX'
To delete the generated files: 'make clean'
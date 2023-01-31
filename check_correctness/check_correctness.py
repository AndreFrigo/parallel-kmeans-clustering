#Script to check the correctness of the parallel/serial algorithm implemented.
#Run the algorithm implemented and copy the initial centroids in the 'centroids.csv' file in the format given.
#Execute this script from the 'check_correctness' folder by running it giving as parameters the dataset path and k
#The script will return the number of executions (use it to stop the algorithm implemented at the same point)
#Run the algorithm implemented and compare the results obtained
import csv
import re
import math
import sys
from sklearn.cluster import KMeans



ncluster = None
inputFile = None
dim = None


if(len(sys.argv) == 3):
    inputFile = str(sys.argv[1])
    try:
        ncluster = int(sys.argv[2])
    except:
        sys.exit("Wrong arguments!\nGive the following: inputFile, k")
else:
    sys.exit("Wrong number of arguments!\nGive the following: inputFile, k")


data = []
with open(inputFile, 'r') as file:
    reader = csv.reader(file)
    cont = 0
    for row in reader:
        # print(row)
        if cont > 0: data.append(tuple(float(e) for e in row))
        cont +=1
    # print (data)
    dim = len(data[0])

centroids = []
c = []
with open('centroids.csv', 'r') as file:
    reader = csv.reader(file)
    cont = 0
    x = None
    v = []
    for row in reader:
        x = row
        x = re.sub('[ ]+', ',', x[0])
        a = x.split(',')
        a = [float(e) for e in a]
        for e in a:
            v.append(e)
    
    par = []
    for i in range(0, len(v)+1):
        if i%dim==0 and i>0:
            c.append(par)
            par = []
        if i<len(v):
            par.append(v[i])
print("CENTROIDS PARALLEL:")
for elem in c:
    print(elem)
print("\n--------------------------------------------------\n")

kmeans = KMeans(n_clusters=ncluster, init=c, n_init=1)
kmeans.fit(data)   
labels = kmeans.predict(data)
centroids  = kmeans.cluster_centers_  
print("NITER: "+str(kmeans.n_iter_)+"\n")

print("CENTROIDS SKLEARN")
for e in centroids: print(e)

import csv
import sys
import random

if len(sys.argv)<3:
    sys.exit("Wrong number of arguments passed, pass number of column and number of entries")

numcol = int(sys.argv[1])
numelem = int(sys.argv[2])

header = []
for i in range(0, numcol):
    header.append("ATTR"+str(i+1))

data = []
i=0
maxran = numelem/10
while(i<numelem):
    el = [round(random.uniform(0, maxran), 2) for e in range(0, numcol)]
    if el not in data:
        data.append(el)
        i+=1


with open('dataset_'+str(numcol)+'_'+str(numelem)+'.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    writer.writerow(header)

    # write the data
    for elem in data:
        writer.writerow(elem)
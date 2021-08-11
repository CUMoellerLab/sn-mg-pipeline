#!/usr/bin/env python

import numpy as np
import sys
from skbio.stats.distance import DistanceMatrix
from prototypeSelection import prototype_selection_destructive_maxdist

print("Running the run_prototypeSelection.py script.")

# import sourmash distance matrix and labels
df = np.load(snakemake.input[0])
labels = open(snakemake.input[1]).read().split()
df2 = 1 - df # similarity to dissimilarity matrix
print("The imported distance matrix has " + str(df2.shape[1]) + " elements.")
pt_min = 2
pt_max = len(labels) - 1

print("Selecting " + str(pt_min) + " to " + str(pt_max) + " prototypes.")

# convert numpy array to type distance matrix
dm = DistanceMatrix(df2)

with open(snakemake.output[0], 'w') as f:
    for k in range(pt_min, pt_max + 1):

        # run prototypeSelection function
        prototypes=prototype_selection_destructive_maxdist(dm, k)
        prototype_labels=", ".join(labels[int(x)] for x in prototypes)

        print(str(k) + "\t" + prototype_labels, file=f)
f.close()

with open(snakemake.log[0], "w") as logfile:
    logfile.write("Running the run_prototypeSelection.py script.\n")
    logfile.write("The imported distance matrix has " + str(df2.shape[1]) + " elements.\n")
    logfile.write("Selecting " + str(pt_min) + " to " + str(pt_max) + " prototypes.\n")
logfile.close()

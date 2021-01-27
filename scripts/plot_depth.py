#!/usr/bin/python
#-------------------
# Author: Yan Hui
# huiyan@food.ku.dk
# 27/01/2021
#-------------------

import os
import pandas as pd
import matplotlib.pyplot as plt

# set the path to the folder containing the .depth files
PATH = '/mnt/md0/RAW_DATA/METAGENOME/results/depth/NXT055/' 

# fetch all files in path
fileNames = os.listdir(PATH)

# filter file name list for files ending with .depth
fileNames = [file for file in fileNames if '.depth' in file]

# loop over all files
for file in fileNames:
    # record the sample name
    sample_name = file.rstrip("\.depth") 
    # read file and append to list
    df = pd.read_csv(PATH + file, sep='\t', names=["genome","position", sample_name], header=None)
    # create line for every file
    plt.plot('position', sample_name, data=df)

# generate the plot
# Put a legend to the right of the current axis
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()
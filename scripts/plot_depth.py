#!/usr/bin/python
#-------------------
# Author: Yan Hui
# E-mail: huiyan@food.ku.dk
# Date: 27/01/2021
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

def normalize_avr_cov(fileName):
    # read file and append to list
    df = pd.read_csv(PATH + fileName, sep='\t', names=['genome','position', 'depth'], header=None)
    # normalize depth by average coverage
    df['normalized'] = df['depth']/df['depth'].mean()
    return df[['position','normalized']]

# calculate the normalized depth of reference
REF = '/mnt/md0/RAW_DATA/METAGENOME/results/depth/NXT055/XIC16.depth'
refName = os.path.basename(REF)
df_ref = normalize_avr_cov(fileName=refName)
df_ref.to_csv('/mnt/md0/RAW_DATA/METAGENOME/results/depth/NXT055/normalized_avr_cov/' + refName, sep = '\t')

# loop over all files
fig, ax = plt.subplots()
for file in fileNames:
    if file == refName:
        continue 
    df = normalize_avr_cov(fileName=file)
    df.to_csv('/mnt/md0/RAW_DATA/METAGENOME/results/depth/NXT055/normalized_avr_cov/' + file, sep = '\t')
    # normalize depth by ref genomes
    df = df.set_index(['position']).div(df_ref.set_index('position')).reset_index()
    
    # rename by sample
    sample_name = file.rstrip('\.depth')
    df.columns = ['position', sample_name]
    df.to_csv('/mnt/md0/RAW_DATA/METAGENOME/results/depth/NXT055/normalized_by_ref/' + file, sep = '\t')
    ax.plot('position', sample_name, data=df)

# generate the plot
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.show()

fig1, ax1 = plt.subplots()
for file in fileNames:
    #if file == refName:
    #    continue 
    df = normalize_avr_cov(fileName=file)
    # normalize depth by ref genomes
    #df = df.set_index(['position']).div(df_ref.set_index('position')).reset_index()
    
    # rename by sample
    sample_name = file.rstrip('\.depth')
    df.columns = ['position', sample_name]
    ax1.plot('position', sample_name, data=df)

# generate the plot
# Put a legend to the right of the current axis
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()
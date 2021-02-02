#!/usr/bin/python
# The python script normalize the base depth by average base depth.
# It also calculates the average depth by sliding windows 
# The results of raw depth and normalzied depth are visualized in the a multi-line chart.
#-------------------
# Author: Yan Hui
# E-mail: huiyan@food.ku.dk
# Date: 27/01/2021
#-------------------
__doc__ = "some function for down-stream processing the depth files"

import matplotlib
matplotlib.use('Agg')
import os
import shutil
import pandas as pd
import matplotlib.pyplot as plt
from itertools import islice
from numpy import mean
import numpy as np
import argparse

def parse_arguments():
    """Read arguments from the console"""
    parser = argparse.ArgumentParser(description="Note: Plot depth files.")
    parser.add_argument("-d", "--depth_dir", help='depth file directory')
    parser.add_argument("-r", "--result_dir", help='output directory')
    parser.add_argument("-f", "--ref", help='one depth file as ref')
    parser.add_argument("-w", "--window_size", help='sliding window size')

    args = parser.parse_args()
    return args

def window(seq, n=2):
    """Returns a sliding window (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...
    """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

# normalize by average
def normalize_avr_cov(depth_dir, fileName):
    df = pd.read_csv(depth_dir + '/' + fileName, sep='\t', names=['genome','position', 'depth'], header=None)
    df['normalized'] = df['depth']/df['depth'].mean()
    return df[['position','normalized']]

# normlazied depth by sliding window
def normalize_sw(seq, window_size):
    rows = []
    i = 1
    for w in window(seq, window_size):
        rows.append([i, sum(w)/len(w)])
        i += 1
    return pd.DataFrame(rows, columns=['windows', 'average depth'])

# plot the raw depth
def raw_depth(fileNames, depth_dir, result_dir):
    fig, ax = plt.subplots()
    # loop over all files
    for file in fileNames:
        # rename by sample
        sample_name = file.rstrip('\.depth')
        
        # read file and append to list
        df = pd.read_csv(depth_dir + '/' + file, sep='\t', names=['genome','position', sample_name], header=None)
        ax.plot('position', sample_name, data=df)
    
    colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired   
    colors = [colormap(i) for i in np.linspace(0, 0.95,len(ax.lines))]
    for i,j in enumerate(ax.lines):
        j.set_color(colors[i])
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # generate the plot
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.savefig(result_dir + '/raw_depth.pdf')

# plot the normalized depth (by average depth)
def normalize_depth(fileNames, depth_dir, result_dir):
    tsv_dir = result_dir + '/normalized_avr_cov/'
    if os.path.exists(tsv_dir):
        shutil.rmtree(tsv_dir)
    os.makedirs(tsv_dir)
    
    fig, ax = plt.subplots()
    # loop over all files
    for file in fileNames:
        df = normalize_avr_cov(depth_dir=depth_dir, fileName=file)
        df.to_csv(tsv_dir + file, sep = '\t', index=False)
        
        # rename by sample
        sample_name = file.rstrip('\.depth')
        df.columns = ['position', sample_name]
        ax.plot('position', sample_name, data=df)
    
    colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired   
    colors = [colormap(i) for i in np.linspace(0, 0.95,len(ax.lines))]
    for i,j in enumerate(ax.lines):
        j.set_color(colors[i])
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # generate the plot
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.savefig(result_dir + '/normalized_avr_cov.pdf')

def normalize_depth_by_ref(fileNames, depth_dir, result_dir, ref):
    refName = os.path.basename(ref)
    df_ref = normalize_avr_cov(depth_dir=depth_dir, fileName=refName)
    tsv_dir = result_dir + '/normalized_by_ref/'
    if os.path.exists(tsv_dir):
        shutil.rmtree(tsv_dir)
    os.makedirs(tsv_dir)
    
    fig, ax = plt.subplots()
    # loop over all files
    for file in fileNames:
        if file == refName:
            continue 
        df = normalize_avr_cov(depth_dir=depth_dir, fileName=file)
        # normalize depth by ref genome
        df = df.set_index(['position']).div(df_ref.set_index('position')).reset_index()
        df.to_csv(tsv_dir + file, sep = '\t', index=False)
        
        # rename by sample
        sample_name = file.rstrip('\.depth')
        df.columns = ['position', sample_name]
        ax.plot('position', sample_name, data=df)
        
    colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired   
    colors = [colormap(i) for i in np.linspace(0, 0.95,len(ax.lines))]
    for i,j in enumerate(ax.lines):
        j.set_color(colors[i])
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.savefig(result_dir + '/normalized_by_ref.pdf')

# sliding windows, window-size=50
def normalized_sliding_windows(fileNames, depth_dir, result_dir, window_size=50):
    window_size = int(window_size)
    tsv_dir = result_dir + '/normalized_Wsize_' + str(window_size) + '/'
    if os.path.exists(tsv_dir):
        shutil.rmtree(tsv_dir)
    os.makedirs(tsv_dir)

    fig, ax = plt.subplots()
    # loop over all files
    for file in fileNames:
        df = normalize_avr_cov(depth_dir=depth_dir, fileName=file)
        df_w = normalize_sw(df['normalized'], window_size)
        df_w.to_csv(tsv_dir + file, sep = '\t', index=False)
        
        # rename by sample
        sample_name = file.rstrip('\.depth')
        df_w.columns = ['windows', sample_name]
        ax.plot('windows', sample_name, data=df_w)
    
    colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired   
    colors = [colormap(i) for i in np.linspace(0, 0.95,len(ax.lines))]
    for i,j in enumerate(ax.lines):
        j.set_color(colors[i])
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # generate the plot
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.savefig(result_dir + '/normalized_Wsize_' + str(window_size) + '.pdf')

def normalized_sw_ref(fileNames, depth_dir, result_dir, window_size, ref):
    window_size = int(window_size)
    
    refName = os.path.basename(ref)
    df_ref = normalize_avr_cov(depth_dir=depth_dir, fileName=refName)
    df_w_ref = normalize_sw(df_ref['normalized'], window_size)
    tsv_dir = result_dir + '/normalized_by_ref_Wsize_' + str(window_size) + '/'
    if os.path.exists(tsv_dir):
        shutil.rmtree(tsv_dir)
    os.makedirs(tsv_dir)

    fig, ax = plt.subplots()
    # loop over all files
    for file in fileNames:
        if file == refName:
            continue 
        df = normalize_avr_cov(depth_dir=depth_dir, fileName=file)
        df_w = normalize_sw(df['normalized'], window_size)
        # normalize depth by ref genome
        df_w = df_w.set_index(['windows']).div(df_w_ref.set_index('windows')).reset_index()
        df_w.to_csv(tsv_dir + file, sep = '\t', index=False)
        
        # rename by sample
        sample_name = file.rstrip('\.depth')
        df_w.columns = ['windows', sample_name]
        ax.plot('windows', sample_name, data=df_w)
    
    colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired   
    colors = [colormap(i) for i in np.linspace(0, 0.95,len(ax.lines))]
    for i,j in enumerate(ax.lines):
        j.set_color(colors[i])
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # generate the plot
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.savefig(result_dir + '/normalized_by_ref_Wsize_' + str(window_size) + '.pdf')

def merge_table(tabs_dir, merged_tab):
    fileNames = os.listdir(tabs_dir)
    fileNames = [file for file in fileNames if '.depth' in file]

    data_merged = []
    df_1 = pd.read_csv(tabs_dir + '/' + fileNames[0], sep='\t')
    data_merged.append(df_1.iloc[:,0])
    sampleNames = []
    for file in fileNames:
        df = pd.read_csv(tabs_dir + '/' + file, sep='\t')
        data_merged.append(df.iloc[:,1])
        sampleName = file.rstrip('\.depth')
        sampleNames.append(sampleName)
    
    colName_1 = list(df_1.columns)
    colNames = sampleNames.insert(0, colName_1[0])
    df_merged = pd.DataFrame(data_merged).T
    df_merged.columns = sampleNames
    df_merged.iloc[:,0] = df_merged.iloc[:,0].astype(int)    
    df_merged.to_csv(merged_tab, sep='\t', index=False)

def merge_table_raw(tabs_dir, merged_tab):
    fileNames = os.listdir(tabs_dir)
    fileNames = [file for file in fileNames if '.depth' in file]

    data_merged = []
    df_1 = pd.read_csv(tabs_dir + '/' + fileNames[0], sep='\t', names=['genome','position', 'depth'], header=None)
    data_merged.append(df_1.iloc[:,1])
    sampleNames = []
    for file in fileNames:
        df = pd.read_csv(tabs_dir + '/' + file, sep='\t', names=['genome','position', 'depth'], header=None)
        data_merged.append(df.iloc[:,2])
        sampleName = file.rstrip('\.depth')
        sampleNames.append(sampleName)
    
    colName_1 = list(df_1.columns)
    colNames = sampleNames.insert(0, colName_1[1])
    df_merged = pd.DataFrame(data_merged).T
    df_merged.columns = sampleNames
    df_merged.to_csv(merged_tab, sep='\t', index=False)

def main():
    args = parse_arguments()
    # fetch all .depth files in path
    fileNames = os.listdir(args.depth_dir)
    fileNames = [file for file in fileNames if '.depth' in file]

    raw_depth(fileNames, args.depth_dir, args.result_dir)
    normalize_depth(fileNames, args.depth_dir, args.result_dir)
    normalize_depth_by_ref(fileNames, args.depth_dir, args.result_dir, args.ref)
    normalized_sliding_windows(fileNames, args.depth_dir, args.result_dir, args.window_size)
    normalized_sw_ref(fileNames, args.depth_dir, args.result_dir, args.window_size, args.ref)

    merge_table(args.result_dir + '/normalized_avr_cov', args.result_dir + '/normalized_avr_cov.tsv')
    merge_table(args.result_dir + '/normalized_by_ref', args.result_dir + '/normalized_by_ref.tsv')
    merge_table(args.result_dir + '/normalized_by_ref_Wsize_' + args.window_size, args.result_dir + '/normalized_by_ref_Wsize_' + args.window_size + '.tsv')
    merge_table(args.result_dir + '/normalized_Wsize_' + args.window_size, args.result_dir + '/normalized_Wsize_' + args.window_size + '.tsv')
    merge_table_raw(args.depth_dir, args.result_dir + '/raw_depth.tsv')

if __name__ == "__main__":
    main()
    
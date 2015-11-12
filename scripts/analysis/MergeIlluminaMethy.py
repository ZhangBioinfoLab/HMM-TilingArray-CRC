# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 23:08:23 2015

@author: Lin
"""

#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-07-05
# email: junjianglin@cs.toronto.com

# This python script is used to merge the illumina DNA methylation rawdata into a 
# standard machine learning matrix
# The data can be reach from NCBI omnibus with accession numberGSE25062
#
# Input List:
# 1. input file list like GSM*.txt
# 2. a pattern used to represent a set of files, like 615453-61557. the number in this range
#    will replace the '*' in input file list

#
# Output List:
# 1. The merged standard machine matrix, rows are samples, and columns are features
#

import pandas as pd
import os

def parseMultiFiles(input, pattern):
    """
    (str,str) => list of str
    Given input files like 'data*.txt', return a list of input files
    """
    index = pattern.split('-')
    input_list = []
    for i in range(int(index[0]),int(index[1])+1):
        input_list.append(input.replace('*',str(i)))
    return input_list
    
def getSpecificDFList(input_list):
    """
    list of str => list of df
    """
    ls_df = []
    for name in input_list:
        df = pd.read_table(name)
        ls_df.append(df.iloc[:,1])
    return ls_df    

def getSpecificDFList_tcga(input_list):
    """
    list of str => list of df
    """
    ls_df = []
    for name in input_list:
        df = pd.read_table(name,skiprows=1)
        ls_df.append(df.iloc[:,1])
    return ls_df   

def main():
    import argparse
    parser = argparse.ArgumentParser(description='summarizing classification statistics')
    parser.add_argument("input",help="input files, must contain * to represent multiple files, or a folder name(all files\
                                      in the folder will be merged")
    parser.add_argument("-p",dest="pattern",help="pattern used for reading multiple\
                        files from input,for example, 1-20 means replacing * with 1 to 20")
    parser.add_argument("-m",dest="mode",default="folder",help="specify the input type from {folder, file},default=[folder]")
    parser.add_argument("output",help="output file prefix")
    
    args = parser.parse_args()
    input = args.input
    output = args.output
    pattern = args.pattern
    mode = args.mode
    
    if mode == "file":
        input_list = parseMultiFiles(input, pattern)
        df_list = getSpecificDFList(input_list)
        df1 = pd.read_table(input_list[0])
        data_matrix = pd.DataFrame(index=df1.iloc[:,0])
        
        for i in range(len(input_list)):
            data_matrix[input_list[i]] = df_list[i].values
        data_matrix = data_matrix.transpose()
        data_matrix.to_csv(output+"_"+pattern+".txt",sep="\t",index=False)
    elif mode == "folder":
        input_list = os.listdir(input)
        input_list = map(lambda x:input+"/"+x,input_list)
        df_list = getSpecificDFList_tcga(input_list)
        df1 = pd.read_table(input_list[0],skiprows=1)
        data_matrix = pd.DataFrame(index=df1.iloc[:,0])
        
        for i in range(len(input_list)):
            data_matrix[input_list[i]] = df_list[i].values
        data_matrix = data_matrix.transpose()
        data_matrix.to_csv(output+".txt",sep="\t",index=False)        
        
        
    
    
    
    
if __name__ == "__main__":
    main()
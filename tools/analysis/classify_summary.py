#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-06-24
# email: junjianglin@cs.toronto.com

# This python script is used to analyze the classification summary data
#
# Input List:
# 1. list of summary files
#
# Output List:
# 1. output summary


import pandas as pd
import numpy as np


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

def getLastModelLine(file_str):
    """ str => int
    Given a file name, return the last model line
    """
    file = open(file_str,'r')
    count = 1
    for line in file:
        if 'model' in line:
            model_count = count
        count = count + 1
    return model_count
    
def getDfList(input_list):
    """ list of str => list of df
    """
    ls_df = []
    for name in input_list:
        ls_df.append(pd.read_table(name,skiprows=getLastModelLine(name)-1,index_col='model'))
    return ls_df
    
def createMultiIndexDF(df_list):
    """ list of df => df
    """    
    m_index = [(x,y) for x in df_list[0].index for y in df_list[0].columns]
    multi_index = pd.MultiIndex.from_tuples(m_index)
    data = []
    for df in df_list:
        data.append(df.values.flatten()[np.newaxis])
    multiIndex_df = pd.DataFrame(data[0],columns=multi_index)
    for i in range(1,len(data)):
        multiIndex_df = multiIndex_df.append(pd.DataFrame(data[i],columns=multi_index),ignore_index=True)
    return multiIndex_df
    
def outputSummary(multiIndex_df,output,pattern):
    """ (multiIndex_df,str,str) => file output
    """
    out_file = open(output+pattern+"_summary.txt","w")
    out_file.write("mean statistics\n")
    out_file.write(str(multiIndex_df.mean())+"\n")
    out_file.write("min statistics\n")
    out_file.write(str(multiIndex_df.min())+"\n")
    out_file.write("max statistics\n")
    out_file.write(str(multiIndex_df.max())+"\n")
    out_file.close()

def plotSummary(multiIndex_df,output,pattern):
    """(multiIndex_df,str,str) => some plot output
    """
    import seaborn as sb
    import matplotlib.pyplot as plt
    new_df = multiIndex_df.swaplevel(0,1,axis=1)
    # drawing boxplot for auc    
    new_df['auc'].boxplot()
    plt.title("Area Under ROC")
    plt.savefig(output+pattern+"_auc.pdf")
    plt.close()
    #drawing boxplot for specificity
    new_df['specificity'].boxplot()
    plt.title("Specificity")
    plt.savefig(output+pattern+"_spec.pdf")
    plt.close()
    #drawing boxplot for sensitivity
    new_df['sensitivity'].boxplot()
    plt.title("Sensitivity")
    plt.savefig(output+pattern+"_sens.pdf")    
    plt.close()
    
    
def main():
    import argparse
    parser = argparse.ArgumentParser(description='summarizing classification statistics')
    parser.add_argument("input",help="input files, must contain * to represent multiple files")
    parser.add_argument("pattern",help="pattern used for reading multiple\
                        files from input,for example, 1-20 means replacing * with 1 to 20")
    parser.add_argument("output",help="output file prefix")
    
    args = parser.parse_args()
    input = args.input
    output = args.output
    pattern = args.pattern
    
    input_list = parseMultiFiles(input, pattern)
    df_list = getDfList(input_list)
    
    ##creating multi-index dataframe    
    multiIndex_df = createMultiIndexDF(df_list)
    ##writing output
    outputSummary(multiIndex_df,output,pattern)
    ##plotting summary
    plotSummary(multiIndex_df,output,pattern)





if __name__ == "__main__":
    main()
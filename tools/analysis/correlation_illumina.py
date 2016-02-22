# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 13:06:47 2015

@author: junjiang
"""

# This python script is used to calculate the pearson coefficient of affymatrix data and illumina methy
# data
#
# Input List:
# 1. Affymatrix data, rows are feature and columns are sample
# 2. illumina merged data, rows are samples and columns are feature
# 3. mapping between Affymatrix and illumina,the format must be:
#    Affy_chr   Affy_start   Affy_end   Affy_name   illu_chr  illu_start   illu_start   illu_end 
#    illu_name  illu***
# 4. output prefix
#
# Output List:
# 1. a scatterplot for the correlation
#


import pandas as pd
import matplotlib as mpl
mpl.use('Agg')    #without using GUI or X server
import seaborn as sb
import matplotlib.pyplot as plt

def cal_cor(cor_df,std = False):
    """
    df => float
    """
    import scipy.stats as st
    
    return st.pearsonr(cor_df.iloc[:,0],cor_df.iloc[:,1])

    
    
def plot_cor(cor_df,std = False):

    zscore = lambda x: (x - x.mean()) / x.std()
    if std:
        cor_df_std = zscore(cor_df)
        cor_df_std.plot(kind="scatter",x='affy',y='illu')
    else:
        cor_df.plot(kind="scatter",x='affy',y='illu')
    



def main():
      import argparse
      parser = argparse.ArgumentParser(description='calculate correlation between Affy and Illu')
	
      parser.add_argument("input1",help="first input file,the result or raw data from Affy")
      parser.add_argument("input2",help="second input file,the result or raw data from illu")
      parser.add_argument("mapping",help="mapping file between Affy and Illu,follow a specific format(see top of this script)")
      parser.add_argument("output",help="output prefix for report and scatterplot")
      args = parser.parse_args()
      
      input1 = args.input1
      input2 = args.input2
      mapping = args.mapping
      output = args.output
      
      
      ## read and preprocessing files
      affy = pd.read_table(input1)
      affy_mean = affy.mean(axis=1)
      illu = pd.read_table(input2)
      illu_mean = illu.mean(axis=0)
      mapping = pd.read_table(mapping,header=None)
      mapping = mapping[[3,7]]
      mapping.rename(columns={3:'affy_name',7:'illu_name'},inplace=True)
      mapping.sort(columns='affy_name',inplace=True)
      
      ##calculate pearson coefficient
      cor_df = pd.DataFrame({'affy':affy_mean[mapping['affy_name']].values,\
                              'illu':illu_mean[mapping['illu_name']].values})
      cor_df.dropna(axis=0,how='any',inplace=True)
      cor = cal_cor(cor_df)
      plot_cor(cor_df)
      plt.title("pearsonr:"+str(cor))
      plt.savefig(output+".pdf")
      
if __name__ == "__main__":
	main()
      
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 17:15:05 2015

@author: Lin
"""

#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-07-05
# email: junjianglin@cs.toronto.com

# This python script is used to find the overlap between identified marker(bed format) and illumina
# Infinium HumanMethylation* DNA methylation assay
#
# Input List:
# 1. input bed file for identified regions
# 2. Infinium HumanMethylation* DNA methylation content file

#
# Output List:
# 1. A list of illumina probe ID(cg*****) which are overlapped with input regions.
#

import pandas as pd






def main():
    import argparse
    parser = argparse.ArgumentParser(description='looking for the overlap between input1 and input2')
    parser.add_argument("input1",help="input1 file, must be in bed format")
    parser.add_argument("input2",help="input2 file, must be in xlsx format, and has the same format as illumina methylation Content file")
    parser.add_argument("output",help="output file prefix")
    
    args = parser.parse_args()
    input1 = args.input1
    input2 = args.input2
    output = args.output
    
    ##### Parsing input file ####
    bed = pd.read_table(input1,header=None,names=['chr','start','end'])
    bed['chr'] = bed['chr'].apply(lambda x:int(x[3:]))
    illumina = pd.read_excel(input2,0)
    
    ### go through illumina and calculate overlap ###
    target_list = []
    for row in illumina.iterrows():
        if row[1]['Chr'] in [4,15,18,20]:
            for row2 in bed.iterrows():
                if row[1]['Chr'] == row2[1]['chr'] and row[1]['MapInfo'] >= row2[1]['start'] and row[1]['MapInfo'] <= row2[1]['end']:
                    target_list.append(row[1]['Name'])
    target = pd.DataFrame(target_list)
    target.to_csv(output+"_overlapIlluminaMethy27.txt",index=False,header=None)
    
    
    
if __name__ == "__main__":
    main()

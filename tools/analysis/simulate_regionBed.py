# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 21:50:52 2015

@author: Lin
"""

#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-06-24
# email: junjianglin@cs.toronto.com

# This python script is used to simulate the identified regions 
#
# Input List:
# 1. input bed file
# 2. number of simulation
# 3. output prefix


#
# Output List:
# 1. simulated bed files
#

import pandas as pd
import random


def write_sim_line(out_file,line_split):
    """ (file,list of str) => write line to file
    """
    # simulation background #chr4:  6137-191260274,   chr15:18303913-100336886,  
    #                        chr18:   1209-76116130,   chr20:8657-62432163
    chr = line_split[0]
    start = int(line_split[1])
    end = int(line_split[2])
    out_file.write(chr +"\t")
    if chr == 'chr4':
        rand_start = random.randint(6137,191260274)
        rand_end = rand_start + (end - start)
        out_file.write(str(rand_start)+"\t"+str(rand_end)+"\n")
    elif chr == 'chr15':
        rand_start = random.randint(18303913,100336886)
        rand_end = rand_start + (end - start)
        out_file.write(str(rand_start)+"\t"+str(rand_end)+"\n")       
    elif chr == 'chr18':
        rand_start = random.randint(1209,76116130)
        rand_end = rand_start + (end - start)
        out_file.write(str(rand_start)+"\t"+str(rand_end)+"\n")
    elif chr == 'chr20':
        rand_start = random.randint(8657,62432163)
        rand_end = rand_start + (end - start)
        out_file.write(str(rand_start)+"\t"+str(rand_end)+"\n")
        
def main():
    import argparse
    parser = argparse.ArgumentParser(description='summarizing classification statistics')
    parser.add_argument("input",help="input files, must contain * to represent multiple files")
    parser.add_argument("-n",dest="number",help="the number of simulation files",type=int,default=1)
    parser.add_argument("output",help="output file prefix")
    
    args = parser.parse_args()
    input = args.input
    number = args.number
    output = args.output
    
    file = open(input,'r')
    for i in range(number):
         random.seed(i+1)
         print "generating {0}th simulated file".format(str(i+1))
         out_file = open(output+"sim"+str(i+1)+".bed",'w')
         for line in file:
             line_split = line.split()
             write_sim_line(out_file,line_split)
         out_file.close()
         file.seek(0)
             
        
if __name__ == "__main__":
    main()

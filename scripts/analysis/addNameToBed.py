# -*- coding: utf-8 -*-
"""
Created on Wed Aug 5 2015

@author: junjiang Lin
@contact: junjianglin@cs.toronto.edu
"""

# This python script is used to add a name column to bed file in the format 'chr,start-end'
#
# Input List:
# 0. bed file
# 1. output prefix

#
# Output List:
# 1. a bed file having a column like chr,start-end
#


import pandas as pd

def main():
	import argparse
	parser = argparse.ArgumentParser(description='generate bed file with chr,start-end')
	parser.add_argument("regions",help="regions to annotate,bed format")
	parser.add_argument("output",help="prefix for output file")
	args = parser.parse_args()
	regions = args.regions
	output = args.output

	bed_df = pd.read_table(regions,header=None)
	bed_df[3] = bed_df[0] + ',' + map(str,bed_df[1]) + '-' +map(str,bed_df[2])
	bed_df = bed_df.iloc[:,:4]
	bed_df.to_csv(output,sep="\t",header=None,index=None)

if __name__ == "__main__":
	main()
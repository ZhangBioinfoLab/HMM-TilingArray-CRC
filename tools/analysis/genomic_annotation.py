# -*- coding: utf-8 -*-
"""
Created on Wed Aug 5 2015

@author: junjiang Lin
@contact: junjianglin@cs.toronto.edu
"""

# This python script is used to summarize the genomic annotation
#
# Input List:
# 0. regions to annotate, bed format
# 1. exon mapping file, -wa -wb in bedtools
# 2. intron mapping file, -wa -wb in bedtools
# 3. promotor mapping file, -wa -wb in bedtools
# 4. intergenic mapping file, -v in bedtools
# 5. output prefix
#
# Output List:
# 1. a tab-delimited file summarize all information
#


import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')    #without using GUI or X server
import seaborn as sb
import matplotlib.pyplot as plt


def main():
	import argparse
	parser = argparse.ArgumentParser(description='summarize mapping files')
	parser.add_argument("regions",help="regions to annotate")
	parser.add_argument("-e","--exon",help="the mapping file for exon generated from bedtools")
	parser.add_argument("-i","--intron",help="the mapping file for intron generated from bedtools")
	parser.add_argument("-g","--intergenic",help="the mapping file for intergenic generated from bedtools")
	parser.add_argument("-p","--promotor",help="the mapping file for promotor generated from bedtools")
	parser.add_argument("output",help="prefix for output file")

	args = parser.parse_args()
	regions = args.regions
	exon = args.exon
	intron = args.intron
	intergenic = args.intergenic
	promotor = args.promotor
	output = args.output

	# Start building summary dataframe
	summary_df = pd.read_table(regions,header=None)
	summary_df[3] = summary_df[0] + ',' + map(str,summary_df[1]) + '-' +map(str,summary_df[2])
	summary_df = summary_df.iloc[:,:4]
	summary_df.rename(columns={0:'chr',1:'start',2:'end',3:'name'},inplace=True)
	summary_df.set_index(keys='name',inplace=True)

	# go through each mapping file
	colors = []
	if promotor:
		promotor_df = pd.read_table(promotor,header=None)
		summary_df['promotor'] = promotor_df[3].value_counts()
		colors.append('gold')	
	if exon:
		exon_df = pd.read_table(exon,header=None)
		summary_df['exon'] = exon_df[3].value_counts()
		colors.append('red')
	if intron:
		intron_df = pd.read_table(intron,header=None)
		summary_df['intron'] = intron_df[3].value_counts()
		colors.append('green')
	if intergenic:
		intergenic_df = pd.read_table(intergenic,header=None)
		summary_df['intergenic'] = intergenic_df[3].value_counts()
		colors.append('blue')

	#post processing
	summary_df.fillna(value=0,inplace=True)
	class_dict = dict(zip(range(7),summary_df.columns))
	class_index = np.argmax(np.array(summary_df.iloc[:,3:]),axis=1)+3
	summary_df['class'] =  [class_dict[i] for i in class_index]
	summary_df.to_csv(output+"_annotation.tsv",sep="\t")
	summary_df['class'].value_counts()[['promotor','intron','exon','intergenic']].plot(kind='pie',autopct='%1.1f%%',shadow=True,startangle=45,colors=colors)
	plt.axis('equal')
	plt.title(output+" genomic annotation")
	plt.savefig(output+"_annotation.pdf")



if __name__ == "__main__":
	main()



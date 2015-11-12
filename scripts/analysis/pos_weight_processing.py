#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-06-23
# email: junjianglin@cs.toronto.com

# This python script is used to filter the position-weight file, like gbm_weight file
#
# Input List:
# 1. position-weight file, must be just two columns, the first column must be position, could be like 
#	 'chr4.101088246.101089486', and the second column must be weight
# 2. optional threshold for selecting first n feautres
# 3. optional threshold for selecting features greater than k weight
# 4. optional separator for position column
#
# Output List:
# 1. bed file that satisfy the thresholds
#


import pandas as pd

def processingFeatureName(features,delimiter):
    """(list of str, str) => dataframe
    for example, 'chr4.101088246.101089486' => chr 101088246 101089486
    """
    import re
    chr = []
    start = []
    end = []
    for feat in features:
        feat_split = re.split('['+delimiter+']',feat)
        chr.append(feat_split[0])
        start.append(int(feat_split[1]))
        end.append(int(feat_split[2])+25)
    out_df = pd.DataFrame({'chr':chr,'start':start,'end':end},columns=['chr','start','end'])
    return out_df

def main():
	import argparse
	parser = argparse.ArgumentParser(description='post processing for position-weight file from CRC project')
	
	parser.add_argument("input",help="input position weight file")
        threshold = parser.add_mutually_exclusive_group()
        threshold.add_argument("-n","--number",dest="number",type=int,help="the number of features will be output")
	threshold.add_argument("-w","--weight",dest="weight",type=float,help="the threshold for weight, only features having weight more than this would be retained")
	parser.add_argument("-d","--delimiter",dest="delimiter",default=".",help="delimiter to separate string in position column, DEFAULT ['.']")
        parser.add_argument("-s","--sort",dest="sort",action='store_false',help="whether to sort pos by weight,DEFAULT[FALSE], this assume the file is already sorted")
        parser.add_argument("output",help="output file name prefix")
        args = parser.parse_args()

	input = args.input
	number = args.number
	weight = args.weight
	delimiter = args.delimiter
        sort = args.sort
        output = args.output


	input = pd.read_table(input,header=None,names=['pos','weight'])
      
        if sort == True:
            input.sort(columns='weight',ascending=0,inplace=True)
            
        if number:
            out = input.iloc[0:number,:]
            out_df = processingFeatureName(list(out['pos']),delimiter)
            print "Finish processing and select first {0} features".format(number)
            print "The minimum weight is {0}".format(out.iloc[-1,:]['weight'])
            out_df.to_csv(output+"_n"+str(number)+".bed",sep="\t",header=None,index=False)
        elif weight:
            out = input[input['weight']>weight]
            out_df = processingFeatureName(list(out['pos']),delimiter)
            print "Finish processing and select features having weight greater than {0} ".format(weight)
            print "The number of filtered features is {0}".format(len(out))
            out_df.to_csv(output+"_w"+str(int(weight))+".bed",sep="\t",header=None,index=False)
        else:
            print "convert file to bed format without fitering"
            out_df = processingFeatureName(list(input['pos']),delimiter)
            out_df.to_csv(output+".bed",sep="\t",header=None,index=False)
            
            
      
if __name__ == "__main__":
    main()
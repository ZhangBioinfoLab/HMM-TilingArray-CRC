#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-06-15
# email: junjianglin@cs.toronto.com

# This python script is used to 
#
# Input List:
#
#
#
# Output List:
#
#


import pandas as pd

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

def readMultiFiles(input_list):
    """
    list of str => list of dataframe
    """
    df_list = []
    for input in input_list:
        df_list.append(readFile(input))
    return df_list

def readFile(input):
    """
    str => dataframe
    """
    return pd.read_table(input,index_col=0)
    
def filterQvalue(df,cutoff=0.05):
    """
    (df,float) => list of str
    """
    return df[df.iloc[:,0]<cutoff].index

def commonFeature(df_ls,cutoff=0.05):
    """
    (list of df,float) = > list of str
    """
    goodFeat = []
    for df in df_ls:
        goodFeat.append(filterQvalue(df,cutoff))
    
    commonFeat = goodFeat[0]
    for feat in goodFeat:
        commonFeat = commonFeat & feat
    
    return commonFeat
    
def singleFilePipeline(input,cutoff):
    input_df = readFile(input)
    input_features = filterQvalue(input_df,cutoff)
    return input_features
    
def multiFilePipeline(input,pattern,cutoff):
    """
    (str,str,float)=>
    """
    input_list = parseMultiFiles(input,pattern)
    input_dfs = readMultiFiles(input_list)
    input_common = commonFeature(input_dfs,cutoff)
    return input_common
    
    
def parseProbe2Bed(ls_str,output):
    """
    (list of str, str) => output file
    for example,a list of str like 'Hs:NCBIv36;chr4-53689339'
    """
    import re
    chr = []
    start = []
    end = []
    for feat in ls_str:
        feat_split = re.split('[;-]',feat)
        chr.append(feat_split[1])
        start.append(int(feat_split[2]))
        end.append(int(feat_split[2])+25)
    out_df = pd.DataFrame({'chr':chr,'start':start,'end':end},columns=['chr','start','end'])
    out_df.to_csv(output+".bed",sep="\t",header=False,index=False)

def parseRange2Bed(ls_str,output):
    """
    (list of str, str) => output file
    for example, a list of str like 'chr15,18303913-18318188'
    """
    import re
    chr = []
    start = []
    end = []
    for feat in ls_str:
        feat_split = re.split('[,-]',feat)
        chr.append(feat_split[0])
        start.append(int(feat_split[1]))
        end.append(int(feat_split[2])+25)
    out_df = pd.DataFrame({'chr':chr,'start':start,'end':end},columns=['chr','start','end'])
    out_df.to_csv(output+".bed",sep="\t",header=None,index=False)
    
def parseFeature2Bed(ls_str,output):
    """
    (list of str, str) => output file
    for example,a list of str like 'Hs:NCBIv36;chr4-53689339'
    or 'chr15,18303913-18318188'
    """
    if len(ls_str) == 0:
        print "nothing to output"
        return
    if 'NCBI' in ls_str[0]:
        parseProbe2Bed(ls_str,output)
    else:
        parseRange2Bed(ls_str,output)
        
def main():
	import argparse
	parser = argparse.ArgumentParser(description='post processing for identified features from CRC project')
	
	parser.add_argument("input",help="first input file,or could be a pattern like 'data*.txt',\
						'*'' will be replaced by the value in option -p")
	parser.add_argument("-s","--second",dest="secondInput",help="second input file, or could be a pattern like \
						'data*.txt','*' will be replaced by the value in option -p")
	parser.add_argument("-p","--pattern",dest="pattern",help="pattern used for reading multiple files from first and second input,\
						for example, 1-20 means replacing * with 1 to 20")
        parser.add_argument("output",help="output file prefix")
        
        parser.add_argument("-c","--cutoff",dest='cutoff',default=0.05,type=float,help='\
                            the cutoff for qvalue,default=0.05')
	
	args = parser.parse_args()
      
        input = args.input
        output = args.output
        second = args.secondInput
        pattern = args.pattern
        cutoff = args.cutoff
        
        
        #Determine which workflow
        direction = 0
        if '*' not in input and second == None:
            direction = 1
        elif '*' not in input and second :
            direction = 2
        elif '*' in input and second == None:
            direction = 3
        elif '*' in input and '*' not in second:
            direction = 4
        elif '*' in input and '*' in second:
            direction = 5
          
          
        if direction == 5:
            print "input: two sets of multiple files"
            input_common = multiFilePipeline(input,pattern,cutoff)
            print "finished first set"
            print "the number of common features in first set is {0}".format(len(input_common))
            print "head features:"
            print input_common[:5]
            second_common = multiFilePipeline(second,pattern,cutoff)
            print "finished second set"
            print "the number of common features in second set is {0}".format(len(second_common))
            print "head features:"
            print second_common[:5]
            
            common_feat = input_common & second_common
            ls_common = list(common_feat)
            if len(ls_common) == 0:
                print 'nothing in common'
            else:
                parseFeature2Bed(ls_common,output)
         
        elif direction == 1:
            print "input: one set of single file"
            ls_feat = singleFilePipeline(input,cutoff)
            print "finished first set"
            print "the number of features in first set is {0}".format(len(ls_feat))
            parseFeature2Bed(ls_feat,output)
            
        elif direction == 2:
            print "input: two sets of single file"
            input_feat = singleFilePipeline(input,cutoff)
            print "finished first set"
            print "the number of features in first set is {0}".format(len(input_feat))
            second_feat = singleFilePipeline(second,cutoff)
            print "finished second set"
            print "the number of features in second set is {0}".format(len(second_feat))
            common_feat = input_feat & second_feat
            parseFeature2Bed(common_feat,output)
        
        elif direction == 3:
            print"input: one set of multiple files"
            input_common = multiFilePipeline(input,pattern,cutoff)
            print "finished first set"
            print "the number of features in first set is {0}".format(len(input_common))
            parseFeature2Bed(input_common,output)
            
        elif direction == 4:
            print"input: two set of files, the first is a set of multiple files,the second is \
                  a single file"
            input_common = multiFilePipeline(input,pattern,cutoff)
            print "finished first set"
            print "the number of features in first set is {0}".format(len(input_common))
            second_feat = singleFilePipeline(second,cutoff)
            print "finished second set"
            print "the number of features in second set is {0}".format(len(second_feat))
            common_feat = input_common & second_feat
            parseFeature2Bed(common_feat,output)
        
        
            
            
                
        


if __name__ == "__main__":
	main()





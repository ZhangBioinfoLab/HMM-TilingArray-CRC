#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-03-20
# email: junjianglin@cs.toronto.com

#
# predict dna methylated regions in tiling array for the colorectal project,
# with B-distribution HMM.
#
# Input List:
# 1. Training data: microarray data for training ( Rdata, or txt format)
# 2. Testing data: microarray data 

#
# Output List:
# 1. training data in RData(biobase) and Txt
# 2. testing data  in RData(biobase) and Txt


# This r script is going to 
#------------------------
rm(list=ls(all=TRUE))

#######################   Initialization   ################
library(optparse)

option_list = list(
        make_option(c("-i","--training"),dest = "training",
                    help="the traing set of tilling array(after normalization),RData or txt"),
        make_option(c("-t","--testing"),dest = "testing",
                    help="the testing set of tilling array(after normalization),RData or txt"),
        make_option(c("-o","--output"),dest = "output",
                    help="the primary name of output file"),
        make_option(c("-d","--differentialOutput"),dest = "diffOut",action="store_true",default=F,
                    help="whether output differential analysis based on t test, will have same output format as hmm"),
        make_option(c("-s","--seed"),type="integer",dest = "seed",default=1,
                    help = "the seed for reproducible research,
                    [default] = 1"),
        make_option("--saveMedium",dest = "saveMedium",action="store_true",default=FALSE,
                    help = "save the medium results to analyze hidden states,
                    [default] = FALSE"),
        make_option("--iniWin",type="integer",dest="iniWin",default=15,
                    help = "Initial window size for merging,
                    [default] = 15"),
        make_option("--proportion",type="double",dest="proportion",default=0.8,
                    help = "proportion of positive states in window,
                    [default] = 0.8"),
        make_option(c("-p","--permutation"),dest="permutation",action="store_true",default=FALSE,
                    help = "whether to do permutation test, if set, permutate the order of probes
                    [default] = FALSE"),
        make_option(c("-n","--number"),dest="number",default=2000,
                    help = "the number of features to retain. [DEFAULT]=2000"),
        make_option("--hmm",dest="hmmoption",default="qvalue",
                    help = "the hmm option,choose from qvalue,logFC,tstatisic,[DEFAULT]=qvalue")
        )
opt = parse_args(OptionParser(option_list=option_list))
if(is.null(opt$training) || is.null(opt$output)|| is.null(opt$testing) ){
        print_help(OptionParser(option_list=option_list))
        q(save = "no")
}

#libraries to load
library(HiddenMarkov)
library(tools)
library(markovchain)
library(MASS)
library(Biobase)


cat("\n","Initializing the environment and loading data...","\n")
set.seed(opt$seed)
training = opt$training
testing = opt$testing
output = opt$output
diffOut = opt$diffOut
saveMedium = opt$saveMedium
proportion = opt$proportion
permutation = opt$permutation
iniWin = opt$iniWin
number = opt$number
hmmoption = opt$hmmoption
# load or read training data
if(file_ext(training) == 'txt'){
    train.matrix = as.matrix(read.delim(file = training))
} else if(file_ext(training) == 'RData'){
    tmp = get(load(training))
    if(class(tmp) == "ExpressionSet"){
        train.matrix = exprs(tmp)
    }else{
        train.matrix = tmp
    }
}

# load or read testing data
if(file_ext(testing) == 'txt'){
    test.matrix = as.matrix(read.delim(file = testing))
} else if(file_ext(testing) == 'RData'){
    tmp = get(load(testing))
    if(class(tmp) == "ExpressionSet"){
        test.matrix = exprs(tmp)
    }else{
        test.matrix = tmp
    }
    
}

#################   Doing permutation ###################
if ( permutation == T){
    cat("\n","permutate the data, shuffle the intensity for each sample...","\n")
    numProbe = nrow(train.matrix)
    for( i in 1:ncol(train.matrix)){
            perm = sample(numProbe,numProbe)
            train.matrix[,i] = train.matrix[perm,i]
    }
}


#################   Generating probe info ################
#input: microarry matrix(feature name contains chromosome and position info)
#output: a dataframe contain each probe info, such as chr, pos, pvalue,pvalue_adjust
getProbeInfo = function(raw.matrix=train.matrix){
  name_list = strsplit(rownames(raw.matrix),'[;-]')
  df_name = data.frame(do.call(rbind,name_list))
  trans_raw.matrix = t(raw.matrix)
  crc = subset(trans_raw.matrix,!grepl('CTR',rownames(trans_raw.matrix)))
  ctr = subset(trans_raw.matrix,grepl('CTR',rownames(trans_raw.matrix)))
  probe_len = nrow(raw.matrix)
  p_values_raw = numeric(probe_len)
  log_fc = numeric(probe_len)
  for(i in 1:probe_len){
     p_values_raw[i] = t.test(crc[,i],ctr[,i])$p.value
     log_fc[i] = log2(mean(as.numeric(crc[,i]))/mean(as.numeric(ctr[,i])))
  }

  p_values = p.adjust(p_values_raw,'fdr') 
  probe = data.frame(chr = df_name$X2,
                     pos = as.numeric(as.character(df_name$X3)),
                     pvalue = p_values_raw,
                     pvalue_ad = p_values,
                     log_fc = log_fc)
}



probe = getProbeInfo(train.matrix)
write.table(data.frame(feature=rownames(train.matrix),qvalue = probe$pvalue_ad,log_fc = probe$log_fc),
            file=paste(output,"_pvalue_logfc",".txt",sep=""),quote=F,row.names=F,sep="\t") 
#################   Naive way:  simple, most differential probe fdr p-value<0.05 ##########
saveWelch = function(probe,diffOut,output,number){
    probeOrder_dex = order(probe$pvalue_ad,decreasing=F)
    if(length(probeOrder_dex)<number){
        idxSigProbe = probeOrder_dex
    }else{
        idxSigProbe = probeOrder_dex[1:number]
    }

    if (length(idxSigProbe)>0){
            cat("\n",length(idxSigProbe)," most significant probes were retained","\n")
            trainSet = as.data.frame(t(train.matrix[idxSigProbe,]))             
            testSet = as.data.frame(t(test.matrix[idxSigProbe,]))
            selectedFeature = rbind(strsplit(make.names(colnames(trainSet)),"\\."))
            featureMatrix = as.data.frame(matrix(unlist(selectedFeature),nrow=length(selectedFeature),byrow=T)[,-c(1,2)])
            names(featureMatrix) = c("chromosome","position")
            featureMatrix$idx = idxSigProbe
            write.table(featureMatrix,file=paste(output,"_Diff_Features",".txt",sep=""),quote=F,sep="\t",col.names=T,row.names=F)
                 
            save(trainSet,file=paste(output,"_trainDiff",".RData",sep=""))
            save(testSet,file=paste(output,"_testDiff",".RData",sep=""))
            write.table(trainSet,file=paste(output,"_trainDiff",".txt",sep=""),quote=F,sep="\t")
            write.table(testSet,file=paste(output,"_testDiff",".txt",sep=""),quote=F,sep="\t")
    }
}




#################   Data Analysis  #####################

#--------------------------------------------------------------------------
#split all probes according to chromosome and sequence of probe into sub-sequences 
#whenever the gap between two probes is larger than max.gap
#--------------------------------------------------------------------------
# function:  splitGroup
# param: probe, the probes user want to split
#        max.gap, any gaps more than max.gap will be split
# return: pvalue_list a 3d-list contains all probes pvalue with subsequence

splitGroup = function(probe,max.gap){
        nchr = length(levels(probe$chr))
        pos_groupByChr = with(probe,split(pos,chr))
        pvalue_groupByChr = with(probe,split(pvalue,chr))
        pvalue_list = lapply(1:nchr,function(x) list())
        names(pvalue_list) = names(pos_groupByChr)
        
        for (i in 1:length(pos_groupByChr)){
                pos_gap_logic = abs(diff(pos_groupByChr[[i]])) > max.gap
                gap.idx = which(pos_gap_logic)
                start = c(1,gap.idx+1)
                end = c(gap.idx,length(pos_groupByChr[[i]]))
                tmpList = mapply(function(s,e,data) data[s:e], start, end, 
                                 MoreArgs = list(pvalue_groupByChr[[i]]),SIMPLIFY=F)
                #print(tmpList)
                pvalue_list[i] = list(tmpList)
        }
        pvalue_list
}



#-----------------------------------------
# for each chromosome, we do K-means cluster to split one chromosome into two clusters
# and estimate the two clusters' mean and variance in each chromosome
#-----------------------------------------

kmeans_emission = function(pvalue_groupByChr,nchr){
    pvalue_list_by_cluster = lapply(1:nchr,function(x) list())
    kmeans_result = lapply(1:nchr,function(x) list())
    Pi_list = lapply(1:nchr,function(x) list())  #transition matrix
    for ( i in 1:nchr ){
            kmeans_result[[i]] = kmeans(pvalue_groupByChr[[i]],centers=2)
            Pi_list[i] = markovchainFit(kmeans_result[[i]]$cluster)
            cluster_index1 = which(kmeans_result[[i]]$cluster == 1)
            cluster_index2 = which(kmeans_result[[i]]$cluster == 2)
            pvalue_list_by_cluster[[i]] = list(pvalue_groupByChr[[i]][cluster_index1],pvalue_groupByChr[[i]][cluster_index2])
    }
    names(pvalue_list_by_cluster) = names(pvalue_groupByChr)
    return(list(emit_list_by_cluster=pvalue_list_by_cluster,Pi_list=Pi_list,kmeans_result=kmeans_result))
}

kmeans_emission_logfc = function(pvalue_groupByChr,nchr){
    pvalue_list_by_cluster = lapply(1:nchr,function(x) list())
    kmeans_result = lapply(1:nchr,function(x) list())
    Pi_list = lapply(1:nchr,function(x) list())  #transition matrix
    for ( i in 1:nchr ){
            kmeans_result[[i]] = kmeans(pvalue_groupByChr[[i]],centers=3)
            Pi_list[i] = markovchainFit(kmeans_result[[i]]$cluster)
            cluster_index1 = which(kmeans_result[[i]]$cluster == 1)
            cluster_index2 = which(kmeans_result[[i]]$cluster == 2)
            cluster_index3 = which(kmeans_result[[i]]$cluster == 3)
            pvalue_list_by_cluster[[i]] = list(pvalue_groupByChr[[i]][cluster_index1],pvalue_groupByChr[[i]][cluster_index2],
                                                pvalue_groupByChr[[i]][cluster_index3])
    }
    names(pvalue_list_by_cluster) = names(pvalue_groupByChr)
    return(list(emit_list_by_cluster=pvalue_list_by_cluster,Pi_list=Pi_list,kmeans_result=kmeans_result))
}





#################  Beta Distribution  #######################

#----------------------------------------------------------------------
# Calculate to get the mean and variance of each cluster in chromosomes
#----------------------------------------------------------------------

#cluster is a data frame that contains all sub cluster information,and it is 
#used to estimate beta parameters
#each row is a chromosome and each column is a indicator of distribution
createQvalueCluster = function(nchr,pvalue_list_by_cluster){
    cluster = as.data.frame(matrix(numeric(nchr*10),nrow = nchr,ncol=10))
    rownames(cluster) = names(pvalue_list_by_cluster)
    colnames(cluster) = c("c1.mean","c1.var","c2.mean","c2.var","c1.size","c2.size",
                          "beta1.alpha","beta1.beta",
                          "beta2.alpha","beta2.beta")
    for ( i in 1:nchr ){
                    cluster[i,]$c1.mean = mean(pvalue_list_by_cluster[[i]][[1]])
                    cluster[i,]$c1.var = var(pvalue_list_by_cluster[[i]][[1]])
                    cluster[i,]$c2.mean = mean(pvalue_list_by_cluster[[i]][[2]])
                    cluster[i,]$c2.var = var(pvalue_list_by_cluster[[i]][[2]])
                    cluster[i,]$c1.size = length(pvalue_list_by_cluster[[i]][[1]])
                    cluster[i,]$c2.size = length(pvalue_list_by_cluster[[i]][[2]])
    }
    return(cluster)
}

createLogFCCluster = function(nchr,logfc_list_by_cluster){
    cluster = as.data.frame(matrix(numeric(nchr*9),nrow = nchr,ncol=9))
    rownames(cluster) = names(logfc_list_by_cluster)
    colnames(cluster) = c("c1.mean","c1.var","c2.mean","c2.var","c3.mean","c3.var","c1.size","c2.size","c3.size")
    for ( i in 1:nchr ){
                    cluster[i,]$c1.mean = mean(logfc_list_by_cluster[[i]][[1]])
                    cluster[i,]$c1.var = var(logfc_list_by_cluster[[i]][[1]])
                    cluster[i,]$c2.mean = mean(logfc_list_by_cluster[[i]][[2]])
                    cluster[i,]$c2.var = var(logfc_list_by_cluster[[i]][[2]])
                    cluster[i,]$c3.mean = mean(logfc_list_by_cluster[[i]][[3]])
                    cluster[i,]$c3.var = var(logfc_list_by_cluster[[i]][[3]])
                    cluster[i,]$c1.size = length(logfc_list_by_cluster[[i]][[1]])
                    cluster[i,]$c2.size = length(logfc_list_by_cluster[[i]][[2]])
                    cluster[i,]$c3.size = length(logfc_list_by_cluster[[i]][[3]])
    }
    return(cluster)
}




#----------------------------------------------------------------------------------
# A function to estimate the parameters of Beta distribution given mean and variance
#-----------------------------------------------------------------------------------
#original way, without using fitdistr function
#estBetaParams <- function(mu, var) {
#        alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
#        beta <- alpha * (1 / mu - 1)
#        return(params = c(alpha = alpha, beta = beta))
#}
#the below change happen in March 22, by Junjiang Lin
estBetaParams <- function(x, mu, var) {
        alpha_start <- ((1 - mu) / var - 1 / mu) * mu ^ 2
        beta_start <- alpha_start * (1 / mu - 1)
        beta_fit = fitdistr(x,"beta",list(shape1=alpha_start,shape2=beta_start)) 
        return(params = c(alpha = beta_fit$estimate[1], beta = beta_fit$estimate[2]))
}

# Step1
#------------------------------------------------
#  Estimate the parameters for beta distribution
#------------------------------------------------

#original way, comment out on Mar. 22
#for ( i in 1:nchr) {
#        estimate1 = estBetaParams(cluster[i,]$c1.mean,cluster[i,]$c1.var)
#        estimate2 = estBetaParams(cluster[i,]$c2.mean,cluster[i,]$c2.var)
#        cluster[i,]$beta1.alpha = estimate1[1]
#        cluster[i,]$beta1.beta = estimate1[2]
#        cluster[i,]$beta2.alpha = estimate2[1]
#        cluster[i,]$beta2.beta = estimate2[2]
#}
estBetaDist = function(pvalue_list_by_cluster,cluster){
    for ( i in 1:nchr) {
        estimate1 = estBetaParams(pvalue_list_by_cluster[[i]][[1]],cluster[i,]$c1.mean,cluster[i,]$c1.var)
        estimate2 = estBetaParams(pvalue_list_by_cluster[[i]][[2]],cluster[i,]$c2.mean,cluster[i,]$c2.var)
        cluster[i,]$beta1.alpha = estimate1[1]
        cluster[i,]$beta1.beta = estimate1[2]
        cluster[i,]$beta2.alpha = estimate2[1]
        cluster[i,]$beta2.beta = estimate2[2]
    }
    return(cluster)
}

estNorDist = function(logFC_list_by_cluster,cluster){

}


#print(cluster)
############################  Hidden Markov Model ########################

#------------------
# hmm function, take in the initial probablities,the cluster, and the index of chromosome i
# output the trained model
#------------------
hmm = function(pvalues,beta1,beta2,Pi,hiddenStates){
        delta = compdelta(Pi)
        x = dthmm(pvalues,Pi,delta=delta,distn="beta",
                        list(shape1 = c(beta1[1],beta2[1]),
                            shape2 = c(beta1[2],beta2[2])))
        #originally without below one line, March 22
        x$y = hiddenStates
        #print(str(x))
        y = BaumWelch(x,bwcontrol(prt=F))
        #y = BaumWelch(x)
        return(y)
}              

hmm_logfc = function(pvalues,normal1,normal2,normal3,Pi,hiddenStates){
        delta = compdelta(Pi)
        x = dthmm(pvalues,Pi,delta=delta,distn="norm",
                        list(mean = c(normal1[1],normal2[1],normal3[1]),
                            sd = c(normal1[2],normal2[2],normal3[2])))
        #x$y = hiddenStates
        print(str(x))
        y = BaumWelch(x)
        return(y)
}

#-----------------------------
#  Basically 3 steps in HMM 
#  1. use the two clusters from each chromosome to estimate initial transition prob(maximum likelihood)
#     and beta parameters (from mean and variance)
#  2. use the estimations from step1 as initial parameter to the whole genome sequence and run EM to get the 
#     trained transition prob and beta parameters
#  3. use trained transition prob and beta para from step2 as initial parameters to subsequences of chromosome, then run EM and 
#     Viterbi to get the predicted states.
#-----------------------------

# Step 2

hmm_train = function(cluster,kmeans_result,pvalue_groupByChr){
    cluster_matrix = as.matrix(cluster)
    hmm_list = lapply(1:nchr,function(x) list())

    for (i in 1:nchr){
        Pi = matrix(nrow=2,ncol=2)
        Pi[1,] = Pi_list[[i]][1]
        Pi[2,] = Pi_list[[i]][2]
        hiddenStates = kmeans_result[[i]]$cluster
        hmm_list[i] = list(hmm(pvalue_groupByChr[[i]],cluster_matrix[i,7:8],
                            cluster_matrix[i,9:10],Pi,hiddenStates))  
    }
    return(list(hmm_list=hmm_list,cluster_matrix=cluster_matrix))
}

hmm_train_logfc = function(cluster,kmeans_result,pvalue_groupByChr){
    cluster_matrix = as.matrix(cluster)
    hmm_list = lapply(1:nchr,function(x) list())

    for (i in 1:nchr){
        Pi = matrix(nrow=3,ncol=3)
        Pi[1,] = Pi_list[[i]][1]
        Pi[2,] = Pi_list[[i]][2]
        Pi[3,] = Pi_list[[i]][3]
        hiddenStates = kmeans_result[[i]]$cluster
        hmm_list[i] = list(hmm_logfc(pvalue_groupByChr[[i]],cluster_matrix[i,1:2],
                            cluster_matrix[i,3:4],cluster_matrix[i,5:6],Pi,hiddenStates))  
    }
    return(list(hmm_list=hmm_list,cluster_matrix=cluster_matrix))
}



# Step 3

hmm_predict = function(nchr,pvalue_list,cluster_matrix){
    pstates_list = lapply(1:nchr,function(x) list())
    chr_names = names(pvalue_list)
    for( i in 1:nchr){
        cat("\n",sprintf("\n%s\n", chr_names[i]))
        posState = ifelse(cluster_matrix[i,1] < cluster_matrix[i,3],1,2)
        backgroundState = ifelse(cluster_matrix[i,1] > cluster_matrix[i,3],1,2)
        Pi = hmm_list[[i]]$Pi
        substates_list = lapply(1:length(pvalue_list[[i]]), function(x) list())
        for ( j in 1:length(pvalue_list[[i]])){
                #if subsequence's length is less than 4, this subsequence will be all the same and determined average pvalue
                if(length(pvalue_list[[i]][[j]]) <= 3){
                        if(mean(pvalue_list[[i]][[j]])<0.2){
                            tmpStates = rep(posState,length(pvalue_list[[i]][[j]]))
                        }else{
                            tmpStates = rep(backgroundState,length(pvalue_list[[i]][[j]]))
                        }
                        
                } else{
                        x = dthmm(pvalue_list[[i]][[j]],Pi,compdelta(Pi),"beta",
                                    pm = hmm_list[[i]]$pm)
                        #y = BaumWelch(x,bwcontrol(prt = F,posdiff = F))
                        y = try(BaumWelch(x,bwcontrol(prt = F)))
                        if(class(y) == "try-error"){
                                tmpStates = Viterbi(x)
                        } else{
                                tmpStates = try(Viterbi(y))
                                if (class(tmpStates) == "try-error"){
                                        tmpStates = Viterbi(x)
                                } 
                        }
                }
                substates_list[j] = list(tmpStates)
        }
        pstates_list[[i]] = c(substates_list,recursive=T)
    }
    return(pstates_list)
}



# update the info in cluster after hmm
update_qvalue_state = function(nchr,pvalue_groupByChr,pstates_list){
    stateInfo = as.data.frame(matrix(numeric(nchr*11),nrow = nchr,ncol=11))
    rownames(stateInfo) = names(pvalue_groupByChr)
    colnames(stateInfo) = c("s1.mean","s1.var","s2.mean","s2.var","s1.size","s2.size",
                          "beta1.alpha","beta1.beta",
                          "beta2.alpha","beta2.beta",
                          "posState")

    pvalue_list_by_state = lapply(1:nchr,function(x) list())
    names(pvalue_list_by_state) = names(pvalue_groupByChr)
    for ( i in 1:nchr ){
            state_index1 = which(pstates_list[[i]] == 1)
            state_index2 = which(pstates_list[[i]] == 2)
            pvalue_list_by_state[[i]] = list(pvalue_groupByChr[[i]][state_index1],pvalue_groupByChr[[i]][state_index2])
            stateInfo[i,]$s1.mean = mean(pvalue_list_by_state[[i]][[1]])
            stateInfo[i,]$s1.var = var(pvalue_list_by_state[[i]][[1]])
            stateInfo[i,]$s2.mean = mean(pvalue_list_by_state[[i]][[2]])
            stateInfo[i,]$s2.var = var(pvalue_list_by_state[[i]][[2]])
            stateInfo[i,]$s1.size = length(pvalue_list_by_state[[i]][[1]])
            stateInfo[i,]$s2.size = length(pvalue_list_by_state[[i]][[2]])
            stateInfo[i,]$beta1.alpha = hmm_list[[i]]$pm$shape1[1]
            stateInfo[i,]$beta1.beta = hmm_list[[i]]$pm$shape2[1]
            stateInfo[i,]$beta2.alpha = hmm_list[[i]]$pm$shape1[2]
            stateInfo[i,]$beta2.beta = hmm_list[[i]]$pm$shape2[2]
            stateInfo[i,]$posState = ifelse(stateInfo[i,]$s1.mean < stateInfo[i,]$s2.mean,1,2)
    }
    return(stateInfo)
}

update_logfc_state = function(nchr,pvalue_groupByChr,pstates_list){
    stateInfo = as.data.frame(matrix(numeric(nchr*7),nrow = nchr,ncol=7))
    rownames(stateInfo) = names(pvalue_groupByChr)
    colnames(stateInfo) = c("s1.mean","s1.var","s2.mean","s2.var","s1.size","s2.size",
                          "posState")

    pvalue_list_by_state = lapply(1:nchr,function(x) list())
    names(pvalue_list_by_state) = names(pvalue_groupByChr)
    for ( i in 1:nchr ){
            state_index1 = which(pstates_list[[i]] == 1)
            state_index2 = which(pstates_list[[i]] == 2)
            pvalue_list_by_state[[i]] = list(pvalue_groupByChr[[i]][state_index1],pvalue_groupByChr[[i]][state_index2])
            stateInfo[i,]$s1.mean = mean(pvalue_list_by_state[[i]][[1]])
            stateInfo[i,]$s1.var = var(pvalue_list_by_state[[i]][[1]])
            stateInfo[i,]$s2.mean = mean(pvalue_list_by_state[[i]][[2]])
            stateInfo[i,]$s2.var = var(pvalue_list_by_state[[i]][[2]])
            stateInfo[i,]$s1.size = length(pvalue_list_by_state[[i]][[1]])
            stateInfo[i,]$s2.size = length(pvalue_list_by_state[[i]][[2]])
            stateInfo[i,]$posState = ifelse(stateInfo[i,]$s1.mean > stateInfo[i,]$s2.mean,1,2)
    }
    return(stateInfo)    
}




# visualize hidden states distribution
# statePos_idx = mapply(function(x,posState) which(x == posState), pstates_list,stateInfo$posState)
# diff_statePos_idx = sapply(statePos_idx,diff)
# diff_all = c(diff_statePos_idx,recursive=T)
# hist(diff_all)

################################  Merge Stage  ##########################

# Merging strategy:
# we use a sliding window with flexible size, the size will be determined by the proportion 
# of state1 probes within this window. In other words, we will set a minimum proportion 
# of state1 probes,says 80%, and the initial sliding window will increase its size until 
# the proportion of state1 probes goes down to 80%, then merge all probes in this window. 
# If the initial proportion of state1 probes is less than 80%, the window will skip and keep 
# moving.
# step size: 1 probe

mergePosStates = function(states,posState,initWindowSize = 20,stateProportion = 1){
        state_end = length(states)
        mergedStates = list()
        start = 1
        windowSize = initWindowSize
        flag = 0
        while (start + windowSize < state_end){
                if(mean(states[start:(start+windowSize)] == posState) >= stateProportion){
                        windowSize = windowSize + 1
                        flag = 1
                } else{
                        if( flag == 1){
                                mergedStates = c(mergedStates,list(c(start,start+windowSize-1)))
                                start = start+windowSize
                                windowSize = initWindowSize
                                flag = 0
                        } else{
                                start = start + 1
                        }
                }
        }
        #whether start to the end could be a merged state?
        if(mean(states[start:state_end] == posState) >= stateProportion){
                mergedStates = c(mergedStates,list(c(start,state_end)))
        }
        return(mergedStates)
}





getMergedDataFrame = function(mergedStates,raw.matrix){
        # Initialize the merged train matrix
        merged_raw.matrix=matrix(NA,nrow=sum(sapply(mergedStates,length)),ncol=ncol(raw.matrix))
        colnames(merged_raw.matrix) = colnames(raw.matrix)
        
        # split raw.matrix by chromosome and get the framework name for merged_train.matrix
        raw.dataframe = data.frame(raw.matrix)
        name_list = strsplit(rownames(raw.matrix),'[;-]')
        df_name = data.frame(do.call(rbind,name_list))
        raw.dataframe$chr =  df_name$X2
        raw.dataframe$pos =  df_name$X3
        raw.dataframe.list = split(raw.dataframe,raw.dataframe$chr)
        merged_name_list = vector("character")
        nsample = ncol(raw.matrix)
        counter = 1
        for( i in 1:nchr ){
                for ( j in 1:length(mergedStates[[i]])){
                        #naming
                        idx1 = mergedStates[[i]][[j]][1]
                        idx2 = mergedStates[[i]][[j]][2]
                        name_start = paste(raw.dataframe.list[[i]][idx1,]$chr,",",raw.dataframe.list[[i]][idx1,]$pos,sep="")
                        name_end = raw.dataframe.list[[i]][idx2,]$pos
                        merged_name = paste(name_start,name_end,sep="-")
                        merged_name_list = c(merged_name_list,merged_name)
                        
                        #calculating mean for merged probes
                        merged_raw.matrix[counter,] = colMeans(as.matrix(raw.dataframe.list[[i]][idx1:idx2,1:nsample]))
                        counter = counter + 1     
                }
        }
        rownames(merged_raw.matrix) = merged_name_list
        merged_raw.dataframe = as.data.frame(merged_raw.matrix)
        return(merged_raw.dataframe)
}




######################### Calculate the merged region's pvalue #############
# after you merge N probes, there will be M regions (M < N). 
# For each region x, you have n_x probes in it. For each sample, take the average 
# signals of the n_x probes to represent signal of the x region. Suppose there are
# T samples, then you will have T signal values for region x. Perform t-test by comparing
# the signal of region x b/w cases and controls. Do this for all M regions.

getSig_merged = function(merged_train.dataframe,number){
    merged_probe_num = nrow(merged_train.dataframe)
    ctr.dataframe = merged_train.dataframe[,grepl("CTR",colnames(merged_train.dataframe))]
    crc.dataframe = merged_train.dataframe[,!grepl("CTR",colnames(merged_train.dataframe))]
    merged_pvalues = numeric()
    merged_logfc = numeric()
    for(i in 1:merged_probe_num){
            merged_pvalues = c(merged_pvalues,t.test(crc.dataframe[i,],ctr.dataframe[i,])$p.value)
            merged_logfc = c(merged_logfc,log2(mean(as.numeric(crc.dataframe[i,]))/mean(as.numeric(ctr.dataframe[i,]))))
    }
    merged_pvalues_ad = p.adjust(merged_pvalues,'fdr') 
    merged_pvalues_order_dec = order(merged_pvalues_ad,decreasing=F)
    if(length(merged_pvalues_order_dec)<number){
        sigIndex = merged_pvalues_order_dec
    }else{
        sigIndex = merged_pvalues_order_dec[1:number]
    }
    sigMergedStates_pvalues_ad = merged_pvalues_ad[sigIndex]
    sigMergedStates_logfc = merged_logfc[sigIndex]
    merged_feature_names = rownames(merged_train.dataframe)
    sigMergedFeatureNames = merged_feature_names[sigIndex]
    return(list(sigIndex=sigIndex,sigMergedStates_pvalues_ad=sigMergedStates_pvalues_ad,sigMergedStates_logfc=sigMergedStates_logfc,
                merged_feature_names=merged_feature_names,sigMergedFeatureNames=sigMergedFeatureNames,merged_pvalues_ad=merged_pvalues_ad,
                merged_logfc = merged_logfc))
}


########################## Generating Training and Testing data for classification ######

getTrain_test_data = function(merged_train.dataframe,sigIndex,mergedStates,output,test.matrix,flag){
    # training set
    trainSet.dataframe = merged_train.dataframe[sigIndex,]
    #transpose the dataframe
    trainSet.dataframe = as.data.frame(t(trainSet.dataframe))

    # testing set
    cat("\n","Start getting the merged data frame of testing set...","\n")
    merged_test.dataframe = getMergedDataFrame(mergedStates,test.matrix)
    testSet.dataframe = merged_test.dataframe[sigIndex,]
    testSet.dataframe = as.data.frame(t(testSet.dataframe))


    cat("\n","Start saving training and testing dataset for machine learning ...","\n")
    trainSet = trainSet.dataframe
    testSet = testSet.dataframe
    selectedFeature = rbind(strsplit(make.names(colnames(trainSet)),"\\."))
    featureMatrix = matrix(unlist(selectedFeature),nrow=length(selectedFeature),byrow=T)
    write.table(featureMatrix,file=paste(output,"_HMM_Features",flag,".txt",sep=""),quote=F,sep="\t",col.names=F,row.names=F)
    save(trainSet,file=paste(output,"_trainHMM",flag,".RData",sep=""))
    save(testSet,file=paste(output,"_testHMM",flag,".RData",sep=""))
    write.table(trainSet,file=paste(output,"_trainHMM",flag,".txt",sep=""),quote=F,sep="\t")
    write.table(testSet,file=paste(output,"_testHMM",flag,".txt",sep=""),quote=F,sep="\t")
}



#########  Driver program ######
cat("\n","generating probe info...","\n")
probe = getProbeInfo(train.matrix)
write.table(data.frame(feature=rownames(train.matrix),qvalue = probe$pvalue_ad,log_fc = probe$log_fc),
            file=paste(output,"_qvalue_logfc",".txt",sep=""),quote=F,row.names=F,sep="\t") 
if(diffOut==T){
    saveWelch(probe,diffOut,output,number)
}
probe_len = nrow(probe)
cat("\n","the size of the data is...",probe_len)
nchr = length(levels(probe$chr))

if(hmmoption == 'qvalue'){
    cat("\n","Start splitting data into subgroups...","\n")
    pvalue_list = splitGroup(probe,1000)
    pvalue_groupByChr = with(probe,split(pvalue,chr))
    #print(str(pvalue_groupByChr))
    #pvalue_groupByChr = with(probe,split(pvalue,chr))
    #print(str(pvalue_groupByChr))

    cat("\n","Start K-mean Clustering...","\n")
    qvalue_kmeans = kmeans_emission(pvalue_groupByChr,nchr)
    pvalue_list_by_cluster = qvalue_kmeans$emit_list_by_cluster
    Pi_list = qvalue_kmeans$Pi_list
    kmeans_result = qvalue_kmeans$kmeans_result
    cluster = createQvalueCluster(nchr,pvalue_list_by_cluster)

    cat("\n","Start estimating beta distribution from the clusters(step1 for hmm)...","\n")
    cluster = estBetaDist(pvalue_list_by_cluster,cluster)
    print(cluster)

    cat("\n","Start estimating chromosome transition prob and beta parameters(step2 for hmm)...","\n")
    qvalue_hmm = hmm_train(cluster,kmeans_result,pvalue_groupByChr)
    hmm_list = qvalue_hmm$hmm_list
    cluster_matrix = qvalue_hmm$cluster_matrix

    cat("\n","Start predicting hidden states (step3 for hmm)...","\n")
    pstates_list = hmm_predict(nchr,pvalue_list,cluster_matrix)

    cat("\n","Start summarizing state info...","\n")
    stateInfo = update_qvalue_state(nchr,pvalue_groupByChr,pstates_list)
    print(stateInfo)

    cat("\n","Start merging positive states...","\n")
    mergedStates = lapply(1:nchr,function(x) list())
    for( i in 1:nchr ){
            mergedStates[i] = list(mergePosStates(pstates_list[[i]],stateInfo$posState[i],
                                                  initWindowSize=iniWin,
                                                  stateProportion=proportion))
    }

    cat("\n","Start getting the merged data frame of training set...","\n")
    merged_train.dataframe = getMergedDataFrame(mergedStates,train.matrix)

    cat("\n","Start calculating the merged region's p value...","\n")
    qvalue_sig_merged = getSig_merged(merged_train.dataframe,number)
    sigIndex = qvalue_sig_merged$sigIndex
    if(saveMedium == T){
            cat("\n","Start saving medium results...","\n")
            merged_pvalues_ad = qvalue_sig_merged$merged_pvalues_ad
            merged_feature_names = qvalue_sig_merged$merged_feature_names
            sigMergedStates_pvalues_ad = qvalue_sig_merged$sigMergedStates_pvalues_ad
            sigMergedStates_logfc = qvalue_sig_merged$sigMergedStates_logfc
            sigMergedFeatureNames = qvalue_sig_merged$sigMergedFeatureNames
            merged_logfc = qvalue_sig_merged$merged_logfc
            write.table(data.frame(feature=merged_feature_names,qvalue = merged_pvalues_ad ,log_fc = merged_logfc),
                file=paste(output,"_hmm_qvalue_logfc",".txt",sep=""),quote=F,row.names=F,sep="\t")
            save(merged_pvalues_ad,merged_feature_names,pstates_list,sigMergedStates_pvalues_ad,sigMergedStates_logfc,sigMergedFeatureNames,
                stateInfo,mergedStates,merged_logfc,file = paste(output,"_mediumResult",".RData",sep=""))
    }

    getTrain_test_data(merged_train.dataframe,sigIndex,mergedStates,output,test.matrix,'qvalue')

}else if(hmmoption == 'logFC'){
    cat("\n","Start splitting data into subgroups...","\n")
    pvalue_list = splitGroup(probe,1000)
    logfc_groupByChr = with(probe,split(log_fc,chr))

    print(str(logfc_groupByChr))

    cat("\n","Start K-mean Clustering...","\n")
    logfc_kmeans = kmeans_emission_logfc(logfc_groupByChr,nchr)
    logfc_list_by_cluster = logfc_kmeans$emit_list_by_cluster
    Pi_list = logfc_kmeans$Pi_list
    kmeans_result = logfc_kmeans$kmeans_result
    cluster = createLogFCCluster(nchr,logfc_list_by_cluster)
    print(cluster)


    cat("\n","Start estimating chromosome transition prob and gaussian parameters(step2 for hmm)...","\n")
    qvalue_hmm = hmm_train_logfc(cluster,kmeans_result,logfc_groupByChr)
    hmm_list = qvalue_hmm$hmm_list
    cluster_matrix = qvalue_hmm$cluster_matrix

    cat("\n","Start predicting hidden states (step3 for hmm)...","\n")
    pstates_list = hmm_predict(nchr,pvalue_list,cluster_matrix)

    cat("\n","Start summarizing state info...","\n")
    stateInfo = update_logfc_state(nchr,pvalue_groupByChr,pstates_list)
    print(stateInfo)


}else if(hmmoption == 'tstatisic'){
    library(tileHMM)
    cat("\n","t statistics: initial estimate...","\n")
    trainLabel = factor(do.call(rbind,strsplit(colnames(train.matrix),'_'))[,2])
    trainLabel_number = sapply(trainLabel,function(x) {if(x=="CRC") return(1) else if(x=="CTR") return(2)})
    train_tstat = shrinkt.st(train.matrix,trainLabel_number,verbose=1)
    hmm_t_init = hmm.setup(train_tstat,state=c("cancer","normal"),pos.state=1)
    gap_split = function(pos_list,max_gap){
         
    }

}



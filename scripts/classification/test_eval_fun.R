# change to your own local path
source('eval_fun.R')

ggout <- 'ROC_test.eps'

data(ROCR.simple)

real.pred <- ROCR.simple$predictions

rand.pred <- ROCR.simple$predictions[sample(length(ROCR.simple$predictions))]

roc.df <- lapply(c("real","rand"), function(x) {
	
	if(x == "real") 
		df <- roceval(real.pred, ROCR.simple$labels, "SVM")
	else if(x == "rand")
		df <- roceval(rand.pred, ROCR.simple$labels, "Random")		
})

roc.df <- do.call(rbind, roc.df)

gg <- evalLinePlot(roc.df, "ROC", "Test ROC")

ggsave(ggout, gg, width=6.4, height=4.8)

system(sprintf("open %s", ggout))


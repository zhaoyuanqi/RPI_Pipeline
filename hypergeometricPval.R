

library(limma)

args <- commandArgs(trailingOnly=T);
argsLen <- length(args);
inputfile <- if (argsLen < 3) 'sample_data.txt' else args[3];


#read the sample data with class annotations (eg control/cancer) at the label column
sampledat<-read.table(inputfile,stringsAsFactors = F)

#generate design matrix
designMat<-model.matrix(~label,data = sampledat)

#fit the linear model
newfit<-lmFit(t(sampledat[,1:ncol(sampledat)-1]),design = designMat)

#Use emperical Bayes statistic for Differential Expression to compute p values, t-statistic, F-statistic, and log odds.
eb_newfit<-eBayes(newfit)

#get list of genes with coeffs and p values
toptable_sample<-topTable(eb_newfit,number = ncol(sampledat)-1)

#set the differentially expressed genes as those with a p value <0.05
diff_genes<-(rownames(toptable_sample[toptable_sample$P.Value<0.05,]))

#total genes for the nypergeometric test, N
total_genes<-rownames(toptable_sample)

#output txt from convertBedToGeneNames.R
boundgenes<-as.matrix(read.table(args[1]))

boundgenes<-boundgenes[,1]

#Only those that are in out final set from the provided experiment
boundgenes_final<-intersect(boundgenes,total_genes)

#calculate the hypergeometric p value
phyper(length(intersect(boundgenes_final,diff_genes)),length(diff_genes),length(total_genes)-length(diff_genes),length(boundgenes_final),lower.tail = F)

#we get a p value of 0.2688612, not significant

#Objective 2 output genes from ensemblToHGNC.R
degenessam68<-as.matrix(read.table(args[2]))

degenessam68<-degenessam68[,1]

degenessam68_final<-intersect(degenessam68,total_genes)

#calculate hypergeometric p value
phyper(length(intersect(degenessam68_final,diff_genes)),length(diff_genes),length(total_genes)-length(diff_genes),length(degenessam68_final),lower.tail = F)

#we get a p value of 0.02416907, significant

outputmat<-matrix(nrow=2,ncol=2)
outputmat[1,1]<-"Pval Intersect of RBP binding genes and experimet diff exp genes : "
outputmat[2,1]<-"Pval Intersect of diff exp genes from RBP knockdown and experiment diff exp genes : "
outputmat[1,2]<-phyper(length(intersect(boundgenes_final,diff_genes)),length(diff_genes),length(total_genes)-length(diff_genes),length(boundgenes_final),lower.tail = F)
outputmat[2,2]<-phyper(length(intersect(degenessam68_final,diff_genes)),length(diff_genes),length(total_genes)-length(diff_genes),length(degenessam68_final),lower.tail = F)
colnames(outputmat)=c("Instance","p.value")

write.table(outputmat,"pvalues.txt",row.names = F,quote = F)
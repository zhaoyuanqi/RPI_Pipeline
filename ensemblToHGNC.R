

library('biomaRt') #needed to convert ensembl ids to gene symbols

args <- commandArgs(trailingOnly=T);

degenes<-as.matrix(read.table(args[1]))
degenes<-degenes[,1]
split_genes<-strsplit(degenes,"\\.")

for (i in 1:length(split_genes)) {
  degenes[i]<-split_genes[[i]][1]
}


mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

mappings <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                          "hgnc_symbol", "description"),values=degenes,mart= mart)
converted_genes<-mappings$hgnc_symbol

for (i in 1:length(converted_genes)) {
  if(converted_genes[i] == ""){
    converted_genes[i] <- NA
  }
}

converted_genes<-na.omit(converted_genes)

write.table(converted_genes,args[2],quote = F, row.names = F,col.names = F)
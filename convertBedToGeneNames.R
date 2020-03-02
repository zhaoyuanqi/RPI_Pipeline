

library(org.Hs.eg.db)  #convert entrez ids to gene symbols

library(CompGO) #Annotate crosslink sites to gene locations using USC Data (txdb) 
                #REQUIRES A JAVA VERSION
                #THAT HAS THE SAME ARCH AS YOUR R VERSION (64 bt R requires 64 bit java)

library(TxDb.Hsapiens.UCSC.hg19.knownGene) #location to gene mappings

library(genomation) #read in bed file as a granges object

args <- commandArgs(trailingOnly=T);
#argsLen <- length(args);
#inputfile <- if (argsLen < 1) 'sample_data.txt' else args[1];




# 'rep1_2_pureclip_crosslink_sites.bed' is the bed file output from pureclip

#load in the bed file as a granges object
grangeBedPureClip<-readBed(args[2], track.line = FALSE, remove.unusual = FALSE,
                           zero.based = TRUE)

grangeBedPureClip_control<-readBed(args[1], track.line = FALSE, remove.unusual = FALSE,
                           zero.based = TRUE)

#transcript database for human genes
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#annotate pureclip file, with a window range of 5kb to cover regions
bedAnnotatedPureClip<-annotateBedFromDb(, gRanges = grangeBedPureClip, db = txdb,
                                        window = 5000)

bedAnnotatedPureClip_control<-annotateBedFromDb(, gRanges = grangeBedPureClip_control, db = txdb,
                                        window = 5000)

entids<-bedAnnotatedPureClip$gene_id@unlistData #entrez ids of annotated bed file

entids_control<-bedAnnotatedPureClip_control$gene_id@unlistData #entrez ids of annotated bed file

idtosymbol<-org.Hs.egSYMBOL
mapped_genes <- mappedkeys(idtosymbol)

#map of entrez ids to gene names
finalMap <- as.list(idtosymbol[mapped_genes])

idstonames<-finalMap[entids]

idstonames_control<-finalMap[entids_control]


geneNamesPureClip<-character(length = 0)
for (i in 1:length(idstonames)) {
  geneNamesPureClip<-base::union(geneNamesPureClip,idstonames[[i]])
}

geneNamesPureClip_control<-character(length = 0)
for (i in 1:length(idstonames_control)) {
  geneNamesPureClip_control<-base::union(geneNamesPureClip_control,idstonames_control[[i]])
}

write.table(base::setdiff(geneNamesPureClip,geneNamesPureClip_control),args[3],quote = F, row.names = F,col.names = F)
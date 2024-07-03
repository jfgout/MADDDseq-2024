################################################################################
#         myDESeq2.R
#
# A set of functions for performing DESeq2 analysis and displaying the 
# results in a Rmarkdown document.



###############################################################################
# jf_runDESeq2
# Performs the DESeq2 differential expression analysis from a pre-loaded txi
# object and a DEseqExperimentSummary "samples" object.
jf_runDESeq2 <- function(samples, tx2gene, design, reduced=NA, MIN_ROW_SUM, cooksCutoff = NA){

  txi <- tximport(samples$file, type="kallisto", tx2gene=tx2gene)
  dds <- DESeqDataSetFromTximport(txi, colData = samples, design = design)
  
  ddsc <- estimateSizeFactors(dds)
  c = counts(ddsc, normalized=TRUE)
  #colnames(c) = paste(samples$mouse_id,samples$genotype,samples$stimulus,sep="_")
  
  keep <- rowSums(counts(dds)) >= MIN_ROW_SUM
  dds <- dds[keep,]
  
  if( is.na(reduced) == T ){
    dds <- DESeq(dds)
  } else {
    dds <- DESeq(dds, test="LRT", reduced = reduced)
  }
  
  res = NA
  if( is.na(cooksCutoff) == T ){
    res = results(dds)
  } else {
    res = results(dds, cooksCutoff = Inf)
  }

  lRes = list(DESeq2res = res, counts = as.data.frame(c), dds=dds )
  lRes
}



################################################################################
# jf_prepareSamples
# This function prepares the data.frame that contains the experiment design 
# information so that it can be used by DESeq2.
# This consists in selecting the time points passed in "vTimes" (useful to remove 
# some time points for some analyzes) and turning all the columns into factors
jf_prepareSamples <- function(designTable, vTimes){
  
  samples = designTable[which( is.element(designTable$Time, vTimes)==T ) ,]
  
  samples$SampleName = substr(samples$FileName, 1, nchar(samples$FileName)-3) # <- This is to remove the ".h5" extension
  samples$file = paste(BD_H5, samples$FileName, sep="/")
  rownames(samples) = paste(samples$SampleName)
  samples$Genotype = factor(samples$Genotype, levels=c("WT", "MGT"))
  samples$Treated = factor(samples$Treated, levels=c("YES", "NO"))
  samples$ReplicationStatus = factor(samples$ReplicationStatus, levels=c("DIVIDING", "ARRESTED"))
  samples$Time = factor(samples$Time, levels=vTimes)
  samples
}




################################################################################
# jf_printDESeqResult

jf_printDESeqResult <- function(lRes, fOut, maxPval=0.05, maxGenesDisplay=20){
  res = lRes[["DESeq2res"]]
  c = lRes[["counts"]]
  dds = lRes[["dds"]]

  sigGenes <- subset(res, padj < maxPval)
  nbSig = nrow(sigGenes)
  nbGenesKept = nrow(c)

  cat("I found a total of ", nbSig, " significant genes with a padj < ", maxPval, "\n", sep="")
  cat("\n  \n  \n")

  to = res[order(res$pvalue, decreasing=F),]
  tor = to[1:maxGenesDisplay , ]
  cat("The top ", maxGenesDisplay, " genes are: \n  \n  \n")
  print(kable(tor))
  
  write.csv(to, file=fOut)
  cat("\n  \n  \nThe full results are available in [", fOut, "](", fOut, ")\n  \n", sep="")

  to
}



#################################################################################
#                                                                               #
#       GENE ONTOLOGY ANALYSIS FUNCTION                                         #
#                                                                               #
#################################################################################


getGeneUniverse <- function(deseqRes){
  vgs = as.numeric(deseqRes$pvalue)
  names(vgs) = rownames(deseqRes)
  vgs
}


####################################################
# my_runGO: a wrapper to run all the topGO analysis
#
# statForGO:  Which statistical method to use in the "runTest" function
# algorithm:  which algorithm to use in the "runTest"
# geneUniverse: a vector containing the values to use for each gene (typically the p-value from DEseq2) with the names of each element being the gene ID
# geneSel:  The function to use for gene selection
# ontology: Which gene ontology to run the analysis on? (One of: BP, MF, CC)
# padjust:  method to use for multiple testing correction. Specify one of: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none" (default is fdr)
# padjustCutOff:  return only categories with an adjust p-value below this cutoff (set to 1 to return all the data)
# MIN_ANNOTATED: remove categories with fewer than MIN_ANNOTATED genes before p-value adjustement
my_runGO <- function( 
  statForGO = c("ks", "fisher", "t", "globaltest", "sum")[1], 
  algorithm = c("classic", "elim", "weight", "weight01", "lea", "parentchild")[1], 
  geneUniverse = c(), 
  geneSel, 
  ontology = "BP", 
  nodeSize=10, 
  annot=annFUN.org, mappingPackage="org.Sc.sgd", ID = "ensembl",
  padjust = "fdr",
  padjustCutOff = 0.05,
  MIN_ANNOTATED = 0
){
  
  
  #resGO = res
  #geneUniverse =  getGeneUniverse(resGO)
  
  allGO2genes <- annFUN.org(whichOnto="ALL", feasibleGenes=NULL, mapping="org.Sc.sgd.db", ID="ensembl")
  GOdata <- new("topGOdata", ontology = ontology, allGenes = geneUniverse , annot = annFUN.GO2genes, 
                GO2genes=allGO2genes, geneSel = sel_pval, nodeSize=5)
  
  #statForGO = "ks"
  #algorithm = "elim"
  #MIN_ANNOTATED = 10
  #padjust = "fdr"
  #padjustCutOff = 0.05
  result = runTest(GOdata, statistic = statForGO, algorithm = algorithm)
  
  allGO = usedGO(object = GOdata)
  allRes <- GenTable(GOdata, pval = result, topNodes = length(allGO))
  allRes = allRes[which(allRes$Annotated>=MIN_ANNOTATED) , ]
  
  allRes$pval = as.numeric(allRes$pval)
  allRes$padj = p.adjust(allRes$pval, method=padjust)
  
  to = allRes[order(allRes$pval, decreasing=F) , ]
  
  ar = allRes[which(allRes$padj<=padjustCutOff) , ]
  
  ar$minv = apply(ar[,c("Significant","Expected")], 1, min)
  ar$maxv = apply(ar[,c("Significant","Expected")], 1, max)
  ar$FE = ar$maxv / ar$minv
  
  ar = ar[ , which( is.element( colnames(ar), c("minv","maxv") ) ==F ) ]
  ar = ar[order(ar$FE, decreasing=T) , ]
  
  for(i in (1:nrow(ar))){
    goID = ar[i,"GO.ID"]
    goTerm = as.character(Term(goID))
    ar[i,"Term"] = goTerm
  }
  
  lres = list(resTable=ar, GOdata=GOdata)
  lres
}
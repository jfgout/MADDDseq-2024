################################################################################
#
# R code to compute mutation rates from MADDD-seq data
#
# The variable BD_DATA should point to a folder with all the output of the 
# MADDD-seq pipeline.
#
# The setwd command should also be modified to point to the required folder on your computer

library(GenomicRanges)
library(Biostrings)
library(BSgenome)
library(genomation)
library(ggplot2)

BD_DATA <- "C:/lab/TR_ERRORS/Projects/YeastMNNG-XP/Jason/data/2022-04-NewPipeLine/yeast"

setwd("C:/lab/TR_ERRORS/Projects/YeastMNNG-XP/Jason/")
source("jason-function-2022-04.R")

vBases = c("A", "T", "C", "G")
v_genotypes = c("WT", "Mgt1")
v_status = c("div", "arr")
v_time_points = c("untx", "t0", "t2", "t6", "t24")

vPrefixStatus = c(Dividing = "div", Arrested = "arr")
vPrefixGenotype = c(WT = "", Mgt1 = "Mgt1-")

MITO_ONLY <- F # <- Set this to true to keep only data from the mitochondrial genome

# Excluding regions with abnormal coverage
tExclude = data.frame(chromosome = c("XII", "XV", "Mito"), start=c(450000,30000,1), end=c(470000,35000,20e6), strand=c("*", "*", "*"))
gExclude = makeGRangesFromDataFrame(tExclude)

fullSpectrum <- data.frame(
  bFrom = character(),
  bTo = character(),
  nbAdduct = numeric(),
  coverage = numeric(),
  adduct = character(),
  rate = numeric(),
  ci_low = numeric(),
  ci_high = numeric(),
  status = character(),
  genotype = character(),
  time_point = character()
)

# To store all the mutations
tmut <- data.frame(chromosome = character(),
                  position = numeric(),
                  from = character(),
                  to = character(),
                  genotype = character(),
                  status = character(),
                  time_point = character()
)


# For debugging purposes:
status <- v_status[1]
genotype <- v_genotypes[2]
time_point <- v_time_points[5]

for(status in v_status){
  for(genotype in v_genotypes){
    for(time_point in v_time_points){
      
      sample_name <- paste(genotype, status, time_point, sep="_")
      cat("Working on ", sample_name, " ...")
      wd <- paste(BD_DATA, "/", genotype, "-", status, "-", time_point, "/filtered_SSC/max_variants_2/", sep="")
      gcov_file <- paste(wd, "coverage.rds", sep = "/")
      #saveRDS(gcov, file = gcov_file)
      gcov <- readRDS(gcov_file)
      
      # !!! WITH THIS NEW VERSION OF THE PIPEPLINE, YOU DO NOT USE THE VCF FILE !!!!
      # !!! ALL MUTATIONS ARE REPORTED IN THE COVERAGE FILE !!!

      
      #mutations_file <- paste(wd, "max_variants_2.DSC.vcf.gz", sep="/")
      #tm <- read.table(mutations_file, h = F)
      #colnames(tm) <- c("chromosome", "start", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "INFO2", "INFO3")
      #tm$chromosome <- substr(tm$chromosome, 4, 10)
      #tm$chromosome[which(tm$chromosome=="M")] <- "Mito"
      #tm$end <- tm$start
      #gad <- makeGRangesFromDataFrame(df= tm, keep.extra.columns = T)
      

      if( MITO_ONLY == T ){
        gcov <- gcov[which(seqnames(gcov)=="Mito")]
        #gad <- gad[which(seqnames(gad)=="Mito")]
      } else {
        gcov <- excludeRegions(gcov, gExclude, ignore.strand = T)
        #gad <- excludeRegions(gad, gExclude, ignore.strand = T)
      }
      
      GRAPH_COVERAGE_AND_MUTAIONS <- F
      if( GRAPH_COVERAGE_AND_MUTAIONS == T ){
        BIN_FACTOR <- 500
        gcov$rpos <- round(start(gcov)/BIN_FACTOR, 0)
        for(chromosome in unique(as.character(seqnames(gcov))) ){
          gcc <- gcov[which(seqnames(gcov)==chromosome)]
          vcov <- by(gcc$depth, gcc$rpos, sum)
          tcov <- data.frame(pos = as.numeric(names(vcov)), cov = as.numeric(vcov))
          tcov$cov <- tcov$cov / BIN_FACTOR
          covFileName <- paste("coverage-", sample_name, "-", chromosome, ".pdf", sep = "")
          pdf(file = covFileName)
          plot(tcov$cov~tcov$pos)
          dev.off(dev.cur())
        }
      }
      
      
      FILTER_SNPS <- F
      if( FILTER_SNPS == T ){
        gcov$rate <- gcov$mutation_count / gcov$depth
        gcov <- gcov[which(gcov$depth>=4 & gcov$rate<3/4)]
        for(base in vBases){
          wm <- which(mcols(gcov)[ , base] > 1)
          mcols(gcov)[wm,base] <- 1
        }
        gcov$mutation_count <- gcov$A + gcov$T + gcov$C + gcov$G
      }
      
      FILTER_FREQUENT_MUTATIONS <- T
      if( FILTER_FREQUENT_MUTATIONS == T ){
        ffreqName <- paste("frequentMutations-", genotype, "_", status, ".rds", sep = "")
        if( MITO_ONLY == T ){
          ffreqName <- paste("frequentMutations-", genotype, "_", status, "-MtONLY.rds", sep = "")
        }
        gu <- readRDS(ffreqName)
        co <- countOverlaps(gcov, gu, ignore.strand = T)
        gcov <- gcov[which(co==0)]
      }
      gcov <- gcov[which(gcov$depth>2 & gcov$mutation_count<2)]
      gmut <- gcov[which(gcov$mutation_count>0)] # This is for the list of mutations that Marc wanted.
      gmut$baseTo <- ""
      for(base in vBases){
        wwb <- which(mcols(gmut)[,base] > 0)
        gmut$baseTo[wwb] <- base
      }
      tmut_tmp <- data.frame(chromosome = as.character(seqnames(gmut)),
                        position = start(gmut),
                        from = gmut$ref,
                        to = gmut$baseTo,
                        genotype = genotype,
                        status = status,
                        time_point = time_point
                        )
      tmut <- rbind(tmut, tmut_tmp)
      
      spectrum <- getSPectrumMutationsRelativeRefStrandFromCovFileOnly(gCoverage = gcov)
      
      spectrum$ci_low <- NA
      spectrum$ci_high <- NA
      for(i in (1:nrow(spectrum))){
        pp <- prop.test(spectrum$nbMut[i], spectrum$coverage[i])
        spectrum$ci_low[i] <- pp$conf.int[1]
        spectrum$ci_high[i] <- pp$conf.int[2]
      }
      spectrum$status <- status
      spectrum$genotype <- genotype
      spectrum$time_point <- time_point
      fullSpectrum <- rbind(fullSpectrum, spectrum)
      cat(" done.\n")
    }
  }
}

write.csv(tmut, file = "all-mutations.csv")

fullSpectrum$time_point <- factor(fullSpectrum$time_point, levels = c("untx", "t0", "t2", "t6", "t24"))




library(dplyr)
myDodge = 1.0
for(genotype in v_genotypes){
  for(status in v_status){
    spec <- fullSpectrum[which(fullSpectrum$genotype == genotype & fullSpectrum$status == status) , ]
    spec <- prepareSpectrumForPlotting(spec)
    spec$time_point <- factor(spec$time_point, levels = c("untx", "t0", "t2", "t6", "t24"))
    myTitle <- paste(genotype, " - ", status, sep = "")
    gp <- ggplot(data = spec, aes(x = mutation, y = rate, fill = time_point)) +
      geom_bar(stat="identity", position=position_dodge()) +
      geom_text(aes(label = nbMut, y = rate), angle = 90, hjust = -0.5, vjust = 0.1, position = position_dodge2(myDodge, preserve = "single")) +
      ggtitle(myTitle)
    plot(gp)
    readline(prompt="Press [enter] to continue")
    
    fOutPlot <- paste("spectrum-mutations-", genotype, "_", status, ".pdf", sep = "")
    if(MITO_ONLY == TRUE){ fOutPlot <- paste("spectrum-mutations-", genotype, "_", status, "-MITO_ONLY.pdf", sep = "") }
    pdf(fOutPlot)
    plot(gp)
    dev.off(dev.cur())
    
    fOutCsv <- paste("spectrum-mutations-", genotype, "_", status, ".csv", sep = "")
    if( MITO_ONLY == TRUE ){ fOutCsv <- paste("spectrum-mutations-", genotype, "_", status, "-MITO_ONLY.csv", sep = "") }
    write.csv(spec, file = fOutCsv)
    
  }
}



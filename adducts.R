################################################################################
#
# R code to compute adduct rates from MADDD-seq data
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
                   strand = character(),
                   from = character(),
                   to = character(),
                   genotype = character(),
                   status = character(),
                   time_point = character()
)



for(status in v_status){
  for(genotype in v_genotypes){
    for(time_point in v_time_points){
      
      sample_name <- paste(genotype, status, time_point, sep="_")
      cat("Working on ", sample_name, " ...")
      wd <- paste(BD_DATA, "/", genotype, "-", status, "-", time_point, "/filtered_SSC/max_variants_2/", sep="")
      cov_file <- paste(wd, "max_variants_2.cov.bin", sep="/")
      #cov_file = "data/Arrested/WT/c1_arr-t0_d1_coverage.bin" # <- Debug to try with old pipeline data
      load(cov_file) # -> loads data into a gcov object
      
      adducts_file <- paste(wd, "max_variants_2.adduct.gtf", sep="/")
      # adducts_file <- "data/Arrested/WT/arr-t0_d1_adducts.gtf" # <- Debug to try with old pipeline data
      gad <- gffToGRanges(gff.file = adducts_file)
      
      tmut_tmp <- data.frame(
        chromosome = as.character(seqnames(gad)),
        position = start(gad),
        strand = as.character(strand(gad)),
        from = gad$adduct,
        to = gad$read_as,
        genotype = genotype,
        status = status,
        time_point = time_point
      )
      tmut <- rbind(tmut, tmut_tmp)
      
      if( MITO_ONLY == T ){
        gcov <- gcov[which(seqnames(gcov)=="Mito")]
        gad <- gad[which(seqnames(gad)=="Mito")]
      } else {
        gcov <- excludeRegions(gcov, gExclude, ignore.strand = T)
        gad <- excludeRegions(gad, gExclude, ignore.strand = T)
      }
      
      
      spectrum <- getSpectrumForAdducts(gc = gcov, adducts = gad)
      
      spectrum$ci_low <- NA
      spectrum$ci_high <- NA
      for(i in (1:nrow(spectrum))){
        pp <- prop.test(spectrum$nbAdduct[i], spectrum$coverage[i])
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

write.csv(tmut, file = "all-adducts.csv", row.names = F)

fullSpectrum$time_point <- factor(fullSpectrum$time_point, levels = c("untx", "t0", "t2", "t6", "t24"))

for(genotype in v_genotypes){
  for(status in v_status){
    spec <- fullSpectrum[which(fullSpectrum$genotype == genotype & fullSpectrum$status == status) , ]
    myTitle <- paste(genotype, " - ", status, sep = "")
    gp <- ggplot(data = spec, aes(x = adduct, y = rate, fill = time_point)) +
      geom_bar(stat="identity", position=position_dodge()) +
      ggtitle(myTitle)
    plot(gp)
    readline(prompt="Press [enter] to continue")
    
    fOutPlot <- paste("spectrum-adducts-", genotype, "_", status, ".pdf", sep = "")
    if(MITO_ONLY == TRUE){ fOutPlot <- paste("spectrum-adducts-", genotype, "_", status, "-MITO_ONLY.pdf", sep = "") }
    pdf(fOutPlot)
    plot(gp)
    dev.off(dev.cur())
    
    fOutCsv <- paste("spectrum-adducts-", genotype, "_", status, ".csv", sep = "")
    if( MITO_ONLY == TRUE ){ fOutCsv <- paste("spectrum-adducts-", genotype, "_", status, "-MITO_ONLY.csv", sep = "") }
    write.csv(spec, file = fOutCsv)
    
  }
}



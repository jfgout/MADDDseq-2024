################################################################################
#
# A collection of useful functions to analyze MADDD-seq data.
#



################################################################################
# add_transcribed_strand
#
# This function adds the information about strand used for transcription to 
# a GRanges object.
add_transcribed_strand <- function(gg, genes_plus, genes_minus) {
  gg$t_strand <- "*"
  co <- countOverlaps(gg, genes_plus, ignore.strand = T)
  gg$t_strand[which(co>0)] = "+"
  co <- countOverlaps(gg, genes_minus, ignore.strand = T)
  gg$t_strand[which(co>0)] = "-"
  gg
}


excludeRegions <- function(gc, gExclude, ignore.strand = T){
  oc = countOverlaps(gc, gExclude, ignore.strand = ignore.strand)
  gc = gc[which(oc==0)]
  gc
}




labelRegions <- function(gg, go){
  gg$isIn = F
  co = countOverlaps(gg, go, ignore.strand=T)
  gg$isIn[which(co>0)] = T
  gg
}




################################################################################
# getSpectrumForAdducts
# Adducts and mutations are treated differently because adducts can happen on either strand!
getSpectrumForAdducts <- function(gc, adducts, vBases = c("A", "T", "C", "G")){
  
  spectrum <- data.frame(
    bFrom = character(),
    bTo = character(),
    nbAdduct = numeric(),
    coverage = numeric()
  )
  #v <- c("A", "T", 2, 10)
  vBasesComp <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  
  i <- 1
  for(bFrom in vBases){
    gcFrom <- gc[which(gc$ref == bFrom | gc$ref == vBasesComp[bFrom])]
    adducstFrom <- adducts[which(adducts$adduct==bFrom)]
    coverage <- sum(gcFrom$depth)
    for(bTo in vBases[which(vBases != bFrom)]){
      nbAdducts <- length(which(adducstFrom$read_as == bTo))
      spectrum[i , ] <- c(bFrom, bTo, nbAdducts, coverage)
      i <- i + 1
    }
  }
  spectrum$nbAdduct <- as.numeric(spectrum$nbAdduct)
  spectrum$coverage <- as.numeric(spectrum$coverage)
  spectrum$adduct <- paste(spectrum$bFrom, spectrum$bTo, sep = "_")
  spectrum$rate <- spectrum$nbAdduct / spectrum$coverage
  spectrum
}



getSPectrumMutationsRelativeRefStrand <- function(gCoverage, gMutations, vBases = c("A", "T", "C", "G")){

  spectrum <- data.frame(
    bFrom = character(),
    bTo = character(),
    nbMut = numeric(),
    coverage = numeric()
  )
  i <- 1
  for(bFrom in vBases){
    gcFrom <- gCoverage[which(gCoverage$ref == bFrom)]
    mutationsFrom <- gMutations[which(gMutations$REF==bFrom)]
    coverage <- sum(gcFrom$depth)
    for(bTo in vBases[which(vBases != bFrom)]){
      nbMutations <- length(which(mutationsFrom$ALT == bTo))
      spectrum[i , ] <- c(bFrom, bTo, nbMutations, coverage)
      i <- i + 1
    }
  }
  spectrum$nbMut <- as.numeric(spectrum$nbMut)
  spectrum$coverage <- as.numeric(spectrum$coverage)
  spectrum$mut <- paste(spectrum$bFrom, spectrum$bTo, sep = "_")
  spectrum$rate <- spectrum$nbMut / spectrum$coverage
  spectrum
  
}

getSPectrumMutationsRelativeRefStrandFromCovFileOnly <- function(gCoverage, vBases = c("A", "T", "C", "G"), collapseMutations = T){

  spectrum <- data.frame(
    bFrom = character(),
    bTo = character(),
    nbMut = numeric(),
    coverage = numeric()
  )
  i <- 1
  totCov <- 0
  totMut <- 0
  for(bFrom in vBases){
    gcFrom <- gCoverage[which(gCoverage$ref == bFrom)]
    coverage <- sum(gcFrom$depth)
    totCov <- totCov + coverage
    for(bTo in vBases[which(vBases != bFrom)]){
      nbMutations <- sum( mcols(gcFrom)[ , bTo])
      if( collapseMutations == T ){
        nbMutations <- length(which(mcols(gcFrom)[,bTo]>0))
      }
      
      totMut <- totMut + nbMutations
      spectrum[i , ] <- c(bFrom, bTo, nbMutations, coverage)
      i <- i + 1
    }
  }
  
  # Adding a row with the total rate:
  spectrum[i,] <- c("N", "N", totMut, totCov)

  spectrum$nbMut <- as.numeric(spectrum$nbMut)
  spectrum$coverage <- as.numeric(spectrum$coverage)
  spectrum$mut <- paste(spectrum$bFrom, spectrum$bTo, sep = "_")
  spectrum$rate <- spectrum$nbMut / spectrum$coverage
  spectrum
  
}



################################################################################
# patternMatchToGRanges
# The function vmatchPattern which I use to find location of triplets does not 
# return a list of GRanges. It returns a list (one entry per scaffold) of data.frame
#
# Parameters:
#   res: The object returned by vmatchPattern
#   strand: one of 2 options ("+" or "-")
#   midPointOnly: Should the GRange only contain the location of the triplet mid point?
patternMatchToGRanges <- function(res, strand, midPointOnly = T){
  
  tres = data.frame(start=numeric(), end=numeric(), width=numeric(), chromosome=character(), strand=character())
  for(seqName in names(res)){
    tmp = as.data.frame(res[[seqName]])
    tmp$chromosome = seqName
    tmp$strand = strand
    tres = rbind(tres, tmp)
  }
  
  if( midPointOnly == T){
    tres$midPoint = tres$start + 1 # !!! CAREFUL: THIS WOULD WORK ONLY FOR TRIPLETS
    tres$start = tres$midPoint
    tres$end = tres$midPoint
    tres$width = 1
  }
  gRes = makeGRangesFromDataFrame(tres)
  gRes
}


################################################################################
# getRatePerTriplet
# Compute the rates for each possible triplet.
# Parameters.
#   adducts: A GRanges object containing the list of adducts to work with.
#             (note: the filtering of G->A should be done before calling this function)
#   gc: A GRanges object containing the information about genome coverage.
#   genome: A BString containing the genome in which triplets should be searched.
#   baseCenter: the focal nucleotide at the center of the triplet (typically G)
getRatePerTriplet_v2 <- function(gam, gc, gTriplet, verbose = F){
  
  co = countOverlaps(gam, gTriplet)
  nbAdducts = length(which(co>0))
  
  # Now looking at coverage
  coc = countOverlaps(gc, gTriplet, ignore.strand=T)
  gcc = gc[which(coc>0)]
  cov = sum(gcc$cov)
  
  vRes = c(nbAdducts, cov)
  vRes
}

################################################################################
# getRateAllTriplets
# This function computes the mutation/adduct rate for each possible triplet (nucleotide context)
# Parameters:
#   gam: a GRanges object with the list of mutations/adducts
#   gc: a GRanges object with covergae information
#   lTriplets: a list of GRanges with the positions of all triplets in the genome (obtained with: getAllTripletsPosition)
#
# Return:
#   A data.frame with mutation/adduct rate per triplet
getRateAllTriplets <- function(gam, gc, lTriplets, focalBase, readAsBase, vBases = c("A", "T", "C", "G")){
  tres <- data.frame(
    focalBase = character(),
    readAs = character(),
    triplet = character(),
    nb = numeric(),
    coverage = numeric()
  )
  iRes <- 0
  for(base5p in vBases){
    for(base3p in vBases){
      triplet = paste(base5p, focalBase, base3p, sep="")
      vres = getRatePerTriplet_v2(gam = gam, gc = gc, gTriplet = lTriplets[[triplet]], verbose = F)
      iRes <- iRes + 1
      tres[iRes,] <- c(focalBase, readAsBase, triplet, vres[1], vres[2])
    }
  }
  tres$nb <- as.numeric(tres$nb)
  tres$coverage <- as.numeric(tres$coverage)
  tres$rate <- tres$nb / tres$coverage
  tres
}


################################################################################
# getAllTripletsPosition
# This function pre-computes the position of each one of the 64 triplets
#
# Parameters:
#   genome  A DNAString object containing the sequence of the genome.
#   vBases: a vector of nucleotides to work with (typically: c("A", "T", "C", "G"))
#
# Return:
#   A list of GRanges (one per triplet)
getAllTripletsPosition <- function(genome, vBases, computeReverseComplement = T, verbose = F){
  lRes = list()
  for(base1 in vBases){
    for(base2 in vBases){
      for(base3 in vBases){
        triplet = paste(base1, base2, base3, sep="")
        
        if( verbose == T){ cat(triplet, "\n") }
        
        DNAtriplet = DNAString(triplet)
        res_f <- vmatchPattern(DNAtriplet, genome)
        gRes_f = patternMatchToGRanges(res_f, "+")
        gRes <- gRes_f
        
        if( computeReverseComplement == T ){
          DNAtriplet_rc = reverseComplement(DNAtriplet)
          res_c <- vmatchPattern(DNAtriplet_rc, genome, fixed=T)
          gRes_c = patternMatchToGRanges(res_c, "-")
          
          gRes = c(gRes_f, gRes_c)
        }
        lRes[[triplet]] = gRes
      }
    }
  }
  lRes
}



# Folds a mutation spectrum
foldSpectrum <- function(spectrum, vBases = c("A", "T", "C", "G"), vbPairs = c(A = "T", T = "A", C = "G", G = "C")){
  
  foldedSpectrum <- data.frame(
    mutation = character(),
    nbMut = numeric(),
    coverage = numeric()
  )
  
  for(bFrom in c("A", "C")){
    for(bTo in vBases[which(vBases != bFrom)]){
      bFromPair <- vbPairs[bFrom]
      bToPair <- vbPairs[bTo]
      foldedMutation <- paste(bFrom, ":", bFromPair , "->", bTo, ":", bToPair, sep = "")
      # cat(foldedMutation, "\n")
      ws <- which( (spectrum$bFrom==bFrom & spectrum$bTo==bTo) | (spectrum$bFrom==bFromPair & spectrum$bTo==bToPair) )
      nbMutations <- sum(spectrum$nbMut[ws])
      coverage <- sum(spectrum$coverage[ws])
      vv <- c(mutation = foldedMutation, nbMut= nbMutations, coverage = coverage)
      foldedSpectrum <- rbind(foldedSpectrum, vv)
    }
  }
  
  w <- which(spectrum$mut == "N_N")
  if( length(w) == 1 ){
    foldedSpectrum <- rbind(foldedSpectrum, c("N->N", spectrum$nbMut[w[1]], spectrum$coverage[w[1]]))
  }
  
  colnames(foldedSpectrum) <- c("mutation", "nbMut", "coverage")
  foldedSpectrum$nbMut <- as.numeric(foldedSpectrum$nbMut)
  foldedSpectrum$coverage <- as.numeric(foldedSpectrum$coverage)
  foldedSpectrum
}

# Adds the 95% confidence intervals to a spectrum
addCI_toFoldedSPectrum <- function(fSpec){
  fSpec$CI_low <- NA
  fSpec$CI_high <- NA
  for(i in (1:nrow(fSpec))){
    nbMut <- fSpec$nbMut[i]
    coverage <- fSpec$coverage[i]
    pp <- prop.test(nbMut, coverage)
    fSpec$CI_low[i] <- as.numeric(pp$conf.int[1])
    fSpec$CI_high[i] <- as.numeric(pp$conf.int[2])
  }
  fSpec
}

# This function calls the spectrum folding function for each time_point in a big spectrum table
prepareSpectrumForPlotting <- function(spectrum){
  ll <- list()
  for(time_point in unique(spectrum$time_point)){
    tSpec <- spectrum[which(spectrum$time_point == time_point) , ]
    foldedSpectrum <- foldSpectrum(spectrum = tSpec)
    foldedSpectrum$time_point <- time_point
    ll[[time_point]] <- foldedSpectrum
  }
  spectrum <- bind_rows(ll)
  spectrum$rate <- spectrum$nbMut / spectrum$coverage
  spectrum <- addCI_toFoldedSPectrum(spectrum)
}

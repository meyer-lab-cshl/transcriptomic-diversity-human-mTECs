# chimSort.R
# March 23 2015 - Modified for LIONS
# ---------------------------------------------------------
# Usage:
#   chimSort.R <ChimericResults> <FilteredOutputName> <Number Exonic Reads>
# 
# Read ChimericReadTool.sh output
#   Parse data into Rdata
#   Calculate derived values
#   Sort/Classify Chimera Cases
#   Output Rdata of sorted chimeric cases
# ---------------------------------------------------------

library(tidyverse)

# CONTROL PANEL ===============================================================

home_directory='/Users/mpeacey/Desktop/thymus-epitope-mapping/ERE-analysis/analysis/External_packages/LIONS/mTEC-analysis'

samples = c('pt214_hi_1', 'pt214_lo_1', 'pt221_hi_1', 'pt221_lo_1', 'pt226_hi_1', 'pt226_lo_1')
lions_csv = list.files(path=home_directory, pattern="*pc.lcsv", full.names=TRUE, recursive=FALSE)
mapped_reads = list.files(path=home_directory, pattern="*mappedReads", full.names=TRUE, recursive=FALSE)

df = data.frame(sample = samples, lions_csv = lions_csv, mapped_reads = mapped_reads) %>%
  mutate(output = glue::glue('{home_directory}/{sample}.lion')) %>%
  mutate(mapped_reads = read_lines(mapped_reads)) %>%
  mutate(scREADS = round(as.numeric(mapped_reads) / 20000000))

for (entry in 1:nrow(df)){
  
  STDIN = c(df[entry, 'lions_csv'],
            df[entry, 'output'],
            df[entry, 'mapped_reads'],
            3, 10, 10, 1, 0.1, 2, 1.5)
  
  INPUT = STDIN[1] # Chimera_Results
  INPUTNAME = unlist(strsplit(as.character(INPUT), split = "\\."))[1]
  INPUTNAME = unlist(strsplit(as.character(INPUTNAME), split = "/"))[11]
  
  OUTPUT = STDIN[2] # Chimeric Output
  
  EXONIC_READS = as.numeric(STDIN[3]) # Number of exonic reads in library
  
  scREADS=max(as.numeric(STDIN[4]), round( EXONIC_READS / 20000000 ) )# >=
  #scREADS=as.numeric(STDIN[4])
  scTHREAD=as.numeric(STDIN[5]) # >=
  scDownThread=as.numeric(STDIN[6]) # >=
  scRPKM=as.numeric(STDIN[7]) # >=
  scCONTR=as.numeric(STDIN[8]) # >=
  scUPCOV=as.numeric(STDIN[9]) # >=
  scUPEXON=as.numeric(STDIN[10])
  
  # IMPORT ======================================================================
  
  # Read output from ChimericReadTool.sh
  ChimeraTable = read.csv(file=INPUT,
                          header=T,
                          sep='\t')
  
  # CALCULATED VALUES ===========================================================
  
  # CONTRIBUTION (via Max)
  # = Repeat(max reads) / Exon(max reads)
  # Relative expression of Repeat to interacting Exon
  
  ChimeraTable$Contribution = (ChimeraTable$RepeatMaxCoverage /
                                 ChimeraTable$ExonMax )
  
  # UpstreamCoverageRatio (UpCov)
  # = Repeat (max) / UpstreamRepeat (max)
  # Look upstream of the repeat to see if there is expression
  
  ChimeraTable$UpCov = (ChimeraTable$RepeatMaxCoverage /
                          ChimeraTable$UpstreamRepeatMaxCoverage)
  
  # Upstream Exon Experssion
  # If Exon 1, set value to Inf
  # = Exon (RPKM) / UpstreamExon (RPKM)
  
  ChimeraTable$UpExonRatio = (ChimeraTable$ExonRPKM /
                                ChimeraTable$UpExonRPKM)
  
  EX1 = which( ChimeraTable[,'exonRankInTranscript'] == 1 )
  ChimeraTable[EX1,'UpExonRatio'] = Inf
  
  # Thread Ratio
  # DownThread / UpThread
  # If UpThread = 0, set value to cutoff
  ChimeraTable$ThreadRatio = (ChimeraTable$DownThread /
                                ChimeraTable$UpThread)
  
  UpThreadZero = which(ChimeraTable$UpThread == 0)
  ChimeraTable[UpThreadZero,'ThreadRatio'] = scTHREAD
  
  # Multi-Exonic Transcripts
  MultiEx = which( ChimeraTable[,'ExonInGene'] == 2)  
  
  # CHIMERA FILTRATION ==========================================================
  
  # Global Filters (Apply to all cases)
  # Initialization
  print('Applying Global Filters to Chimeric Interaction Table')
  print(paste('     Input Entries   : ', nrow(ChimeraTable)))
  
  # CASE 1: Upstream Interactions
  RETAIN_UP = which( ChimeraTable$ER_Interaction == 'Up' &
                       ChimeraTable$Total >= scREADS &
                       ChimeraTable$ThreadRatio >= scTHREAD &
                       ChimeraTable$Contribution >= scCONTR &
                       ChimeraTable$UpCov >= scUPCOV)
  
  # CASE 2: Up Edge Interactions
  RETAIN_UPEDGE = which( ChimeraTable$ER_Interaction == 'UpEdge' &
                           ChimeraTable$Total >= scREADS &
                           ChimeraTable$exonRankInTranscript == 1 &
                           ChimeraTable$ExonInGene == 2 &
                           ChimeraTable$Contribution >= scCONTR/5 &
                           # Order of magnitude reduction since more direct
                           # evidence of promoter is avaiable
                           ChimeraTable$ThreadRatio >= scTHREAD )
  
  # CASE 3: Exon Inside Interactions
  RETAIN_EINSIDE = which( ChimeraTable$ER_Interaction == 'EInside' &
                            ChimeraTable$Total >= scREADS &
                            ChimeraTable$ExonInGene == 2 &
                            ChimeraTable$ThreadRatio >= scTHREAD )
  
  ## CASE 4: Repeat Inside Interaction
  #  RETAIN_RINSIDE = which( ChimeraTable$ER_Interaction == 'RInside' &
  #                          ChimeraTable$Total >= scREADS &
  #			  ChimeraTable$ExonInGene == 2 &
  #			  ChimeraTable$ThreadRatio >= scTHREAD &
  #			  ChimeraTable$DownThread >= scDownThread )
  #
  ## with RInside
  #  ChimeraOut = ChimeraTable[c(RETAIN_UPEDGE, RETAIN_EINSIDE, RETAIN_UP, RETAIN_RINSIDE),]
  
  # without RInside cases
  ChimeraOut = ChimeraTable[c(RETAIN_UPEDGE, RETAIN_EINSIDE, RETAIN_UP),]
  
  # PARSE and OUTPUT ============================================================
  
  print(paste('     Output Entries   : ', nrow(ChimeraOut)))
  # Sort Rows
  SORT=4 # Coordinates
  ChimeraOut = ChimeraOut[order(ChimeraOut[,SORT]),]
  
  # Start of Repeat (For comparisons)
  ChimeraOut$RepeatID = paste('chr',ChimeraOut$Chromosome,':',ChimeraOut$RStart,sep='')
  
  # Add input to last column
  #  (for comparing multiple lists)
  
  #LIBRARY = unlist(strsplit(INPUT,split = '/'))[1]
  LIBRARY = INPUTNAME
  ChimeraOut = cbind(ChimeraOut, LIBRARY)
  
  # Write output CSV
  write.table(ChimeraOut,
              file = OUTPUT,
              quote = F,
              sep = '\t',
              row.names = F,
              col.names = T)
  
}




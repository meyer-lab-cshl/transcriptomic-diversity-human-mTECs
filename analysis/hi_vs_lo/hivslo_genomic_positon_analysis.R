library(GenomicRanges)
library(dplyr)
library(tidyr)

#################################################################
# 
#################################################################

make_GRanges = function(mode, results_df){
  
  if (mode == 'TE'){
    
    annotation = read.table(file = "/Users/mpeacey/TE_thymus/analysis/annotation_tables/hg38_rmsk_TE.gtf.locInd.locations.txt", header = 1)
    annotation = separate(annotation, chromosome.start.stop, into = c('chr', 'start.stop'), sep = ':')
    annotation = separate(annotation, start.stop, into = c('start', 'end'), sep = '-')
    annotation = rename(annotation, locus = TE)
    
    df = merge(results_df, annotation, by = 'locus')
    
  }
  
  if (mode == 'gene'){
    
    annotation = read.table(file = '/Users/mpeacey/TE_thymus/analysis/annotation_tables/gencode.v38_gene_annotation_table.txt', header = 1)
    annotation = select(annotation, c('Geneid', 'Chromosome', 'Start', 'End', 'Strand'))
    annotation = rename(annotation, chr = Chromosome, start = Start, end = End, strand = Strand)
    
    df = merge(results_df, annotation, by = 'Geneid')
    
  }
  
  output = makeGRangesFromDataFrame(df, keep.extra.columns = T)
  
  return(output)
  
}

## GRanges


GRanges_gene = make_GRanges(mode = 'gene',
                            results_df = results_df_local_gene)

GRanges_TE = make_GRanges(mode = 'TE',
                          results_df = results_df_local_TE)



GRanges_gene_sigdiff = make_GRanges(mode = 'gene',
                                    results_df = results_df_local_gene_sigdiff)


GRanges_TE_sigdiff = make_GRanges(mode = 'TE',
                                  results_df = results_df_local_TE_sigdiff)

subsetByOverlaps(GRanges_TE_sigdiff, GRanges_gene_sigdiff)

#################################################################
# regioneR
################################################################### 

random_RS = resampleRegions(GRanges_TE_sigdiff, universe=GRanges_TE)

## Do my differentially expressed TEs overlap with differentially expressed genes more frequently than expected by chance?

pt = permTest(A = GRanges_TE_sigdiff, 
              B = GRanges_gene_sigdiff, 
              ntimes = 100,
              randomize.function = resampleRegions,
              universe = GRanges_TE,
              evaluate.function = numOverlaps)

## Are differentially expressed TEs closer to differentially expressed genes than expected by chance?

pt = permTest(A = GRanges_TE_sigdiff, 
              B = GRanges_gene_sigdiff, 
              ntimes = 100,
              randomize.function = resampleRegions,
              universe = GRanges_TE,
              evaluate.function = meanDistance)

plot(pt) 
library(dplyr)

filename = "/Users/mpeacey/TE_thymus/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
GTEX_data = read.table(filename, sep="\t", nrow=-1, header=TRUE, as.is=TRUE, quote="")

filename = "/Users/mpeacey/TE_thymus/SraRunTable.txt"
SRA_data = read.csv(filename) %>%
  filter(biospecimen_repository == 'GTEx') %>%
  rename(SAMPID = biospecimen_repository_sample_id)

data = merge(GTEX_data, SRA_data, by = 'SAMPID')

testis_data = filter(data, SMTS == 'Testis') %>% 
  filter(LibraryLayout == 'PAIRED') %>%
  filter(SMGEBTCHT == 'TruSeq.v1') %>%
  filter(SMAFRZE == 'RNASEQ') %>%
  filter(SMNABTCHT == 'RNA Extraction from Paxgene-derived Lysate Plate Based') %>%
  filter(SMATSSCR == 0 | 1) %>%
  filter(Is_Tumor == 'No') %>%
  arrange(desc(SMRIN))
testis_data = testis_data[grep("^GCF",testis_data$AssemblyName),]
testis_data = testis_data[1:10, ]

muscle_data = filter(data, SMTS == 'Muscle') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  filter(SMGEBTCHT == 'TruSeq.v1') %>%
  filter(SMAFRZE == 'RNASEQ') %>%
  filter(SMNABTCHT == 'RNA Extraction from Paxgene-derived Lysate Plate Based') %>%
  filter(SMATSSCR == 0 | 1) %>%
  filter(Is_Tumor == 'No') %>%
  arrange(desc(SMRIN))
muscle_data = muscle_data[grep("^GCF",muscle_data$AssemblyName),]
muscle_data = muscle_data[1:10, ]

testis_SRA_IDs = testis_data$Run
muscle_SRA_IDs = muscle_data$Run

write.csv(testis_SRA_IDs, file = '/Users/mpeacey/TE_thymus/testis_data.csv')
write.csv(muscle_SRA_IDs, file = '/Users/mpeacey/TE_thymus/muscle_data.csv')
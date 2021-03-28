library(dplyr)

filename = "/Users/mpeacey/TE_thymus/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
data = read.table(filename, sep="\t", nrow=-1, header=TRUE, as.is=TRUE, quote="")

testis_data = filter(data, SMTS == 'Testis') %>% 
  filter(SMGEBTCHT == 'TruSeq.v1') %>%
  filter(SMAFRZE == 'RNASEQ') %>%
  filter(SMATSSCR == 0 | 1) %>%
  arrange(desc(SMRIN))

muscle_data = filter(data, SMTS == 'Muscle') %>% 
  filter(SMGEBTCHT == 'TruSeq.v1') %>%
  filter(SMAFRZE == 'RNASEQ') %>%
  filter(SMATSSCR == 0 | 1) %>%
  arrange(desc(SMRIN))

filename = "/Users/mpeacey/TE_thymus/SraRunTable.txt"
SRA_data = read.csv(filename)

testis_SRA_IDs = filter(SRA_data, biospecimen_repository_sample_id %in% testis_data[1:10, 1])$Run
muscle_SRA_IDs = filter(SRA_data, biospecimen_repository_sample_id %in% muscle_data[1:10, 1])$Run

write.csv(testis_SRA_IDs, file = '/Users/mpeacey/TE_thymus/testis_data.csv')
write.csv(muscle_SRA_IDs, file = '/Users/mpeacey/TE_thymus/muscle_data.csv')
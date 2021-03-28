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

write.csv(testis_data[1:10, 1], file = '/Users/mpeacey/TE_thymus/testis_data.csv')
write.csv(muscle_data[1:10, 1], file = '/Users/mpeacey/TE_thymus/muscle_data.csv')
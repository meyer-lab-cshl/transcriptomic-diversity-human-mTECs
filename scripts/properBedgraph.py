import sys
import os
import pandas as pd
import re 

#input variables from command line, inputfile and outputfile
bedgraph = sys.argv[1]                  
out_bedgraph = sys.argv[2]              

#convert inputfile into dataframe
df= pd.read_csv(bedgraph, sep='\t', names = ["Chr", "Start", "Coverage", "iterator", "strand", "x"]) 
df= df.drop('x', 1)
df= df.drop('iterator', 1)
df['End'] = df['Start']+1

df['Chr'] = 'chr' + df['Chr'].astype(str)
df.loc[df.Chr == 'chrMT', 'Chr'] = 'chrM'
df = df[df.Chr != 'chrpGIBS-PHE']
df = df[df.Chr != 'chrpGIBS-LYS']
df = df[df.Chr != 'chrpGIBS-THR']

df.loc[df.strand == '-', 'Coverage'] *= -1

df= df.drop('strand', 1)
df = df[['Chr', 'Start', 'End', 'Coverage' ]]

#save dataframe
df.to_csv(out_bedgraph, sep='\t', index=False, header=False)

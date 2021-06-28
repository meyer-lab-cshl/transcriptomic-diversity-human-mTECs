#!/usr/bin/env python

import pandas as pd
import numpy as np

in_dir="/grid/meyer/home/jacarter/TSS/CTSS/Merged_CTSS_BED/"
out_dir="/grid/meyer/home/jacarter/TSS/CTSS/HOMER/"
file="All_thymus_merged.txt"

df_all=pd.read_csv(in_dir+file,sep='\t',header=None,names=['Chr','Start','Stop','Samples','Sample_indices','Strand'])

listing=[]
listing_subject=[]
high=[]
low=[]
for x in df_all['Samples']:
    counter=0
    counter_subject=[]
    no_low=0
    no_hi=0
    for y in x.split(','):
        counter=counter+1
        counter_subject.append(y.split('_')[0])
        if y.split('_')[1]=='lo':
            no_low+=1
        elif y.split('_')[1]=='hi':
            no_hi+=1
    listing.append(counter)
    listing_subject.append(np.unique(counter_subject).shape[0])
    high.append(no_hi)
    low.append(no_low)

df_all['Listing']=np.array(listing)
df_all['Listing_subject']=np.array(listing_subject)
df_all['High']=np.array(high)
df_all['Low']=np.array(low)
df_all['UniqueIdentifier']=range(df_all.shape[0])

in_dir_reads="/grid/meyer/home/jacarter/TSS/CTSS/Paraclu_out/Paraclu_Trim/"
subjects=['pt212_lo', 'pt221_lo', 'pt226_lo','pt87_lo','pt214_lo',
          'pt212_hi', 'pt221_hi', 'pt226_hi','pt87_hi','pt214_hi']
for i,subject in enumerate(subjects):
    x=pd.read_csv(in_dir_reads+subject+'_TPM_paraclu_simplified.txt',sep='\t',header=None)
    x.columns=['Chr','Strand','Start','Stop','Read_Sites','TPM','min','max']
    x=x.assign(idx=x.index+1)
    x=x.assign(sample_name=subject)
    x=x.assign(Name=x.sample_name+'_'+x.idx.astype(str))
    x=x[['TPM','Read_Sites','Name']]
    if i==0:
        df_reads=x
    else:
        df_reads=df_reads.append(x)
    df_all.insert(df_all.shape[1], subject+'_TPM', np.zeros(df_all.shape[0]))
    df_all.insert(df_all.shape[1], subject+'_Read_Sites', np.zeros(df_all.shape[0]))

for i,x in enumerate(df_all.Sample_indices):
    for y in x.split(','):
        tpm,reads=df_reads[df_reads.Name==y][['TPM','Read_Sites']].sum()
        df_all.at[i,y.split('_')[0]+'_'+y.split('_')[1]+'_TPM']=tpm+df_all.iloc[i][y.split('_')[0]+'_'+y.split('_')[1]+'_TPM']
        df_all.at[i,y.split('_')[0]+'_'+y.split('_')[1]+'_Read_Sites']=reads+df_all.iloc[i][y.split('_')[0]+'_'+y.split('_')[1]+'_Read_Sites']

df_all.to_pickle('/grid/meyer/home/jacarter/TSS/df_all.pkl')

threshold=1
df_Low=df_all[(df_all.High==0) & (df_all.Low>(threshold))]
df_Low=df_Low.assign(Class=['Low']*df_Low.shape[0])
df_High=df_all[(df_all.Low==0) & (df_all.High>(threshold))]
df_High=df_High.assign(Class=['High']*df_High.shape[0])

x=df_all[['Chr','Start','Stop','UniqueIdentifier','Samples','Strand']]
x.to_csv(out_dir+'forHomer.txt',sep='\t',header=False,index=False)

x=df_High[['Chr','Start','Stop','UniqueIdentifier','Samples','Strand']]
x.to_csv(out_dir+'Homer_High.txt',sep='\t',header=False,index=False)

x=df_Low[['Chr','Start','Stop','UniqueIdentifier','Samples','Strand']]
x.to_csv(out_dir+'Homer_Low.txt',sep='\t',header=False,index=False)

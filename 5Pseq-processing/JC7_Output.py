#!/usr/bin/env python

import pandas as pd
import numpy as np

df_all=pd.read_pickle('/grid/meyer/home/jacarter/TSS/df_all.pkl')

in_dir_Homer="/grid/meyer/home/jacarter/TSS/CTSS/HOMER/HOMER_Out/"
x=pd.read_csv(in_dir_Homer+'Homerout.txt',sep='\t')
x.index=np.array(x[x.columns[0]])
x=x.sort_values(by=x.columns[0])
x=x.drop(columns=x.columns[0])
features=['Annotation','Distance to TSS','Nearest Ensembl','Gene Name','Gene Description','Gene Type']
df_all[features]=x[features]
df_all['Annotations_short']=[y.split("(")[0] for y in x.Annotation]

x=pd.read_csv(in_dir_Homer+'Homerout_CpG.txt',sep='\t')
x.index=np.array(x[x.columns[0]])
x=x.sort_values(by=x.columns[0])
x=x.drop(columns=x.columns[0])
features=['CpG%','GC%']
df_all[features]=x[features]

threshold=2
df_Low=df_all[(df_all.High==0) & (df_all.Low>=(threshold))]
df_Low=df_Low.assign(Class=['Low']*df_Low.shape[0])
df_High=df_all[(df_all.Low==0) & (df_all.High>=(threshold))]
df_High=df_High.assign(Class=['High']*df_High.shape[0])

df_all.to_pickle('/grid/meyer/home/jacarter/TSS/df_all.pkl')
df_High.to_pickle('/grid/meyer/home/jacarter/TSS/df_High.pkl')
df_Low.to_pickle('/grid/meyer/home/jacarter/TSS/df_Low.pkl')

#get ${HOME}/TSS/df_all*
#get ${HOME}/TSS/CTSS/Power*
#get ${HOME}/TSS/CTSS/HOMER/HOMER_Out/*
#get ${HOME}/TSS/CTSS/HOMER/HOMER_High_Low/*
#get ${HOME}/TSS/CTSS/HOMER/HOMER_Motifs/Motifs_High
#get ${HOME}/TSS/CTSS/HOMER/HOMER_Motifs/Motifs_Low

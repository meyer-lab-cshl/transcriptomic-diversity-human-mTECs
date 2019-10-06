# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Author: jaerveli
# Date: 21 March 2013
# Description:
#      A simple script to collect into one table alignment info from 'uniquelyAligned.info' files
#       and pcr amplification from *_countsPerMol.tsv files.
#
# cat /g/steinmetz/project/mTEC_Seq/src/getPerformace.R | R --slave --args
#       indirAlignment='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/'
#       indirCollapsing='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/'
#	 inAlignmentPattern='uniquelyAligned.info'
#	 inAmplificationPattern='_countsPerMol.tsv'
#       outDir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/performance/'
#       sortingOut='sorting-and-alignment.tsv'
#       amplificationOut='molecules-and-positions.tsv'
#       summaryOut='summary.tsv'
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


options(scipen=20)


#### PROCESS THE INPUT ARGUMENTS

args = commandArgs(TRUE);

if(length(args)==0){
        stop("No arguments provided!\n")
} else {
        indirAlignment =ifelse( length(grep('indirAlignment', args) )>0, strsplit(args[grep('indirAlignment', args)], '=')[[1]][2], stop('No indirAlignment argument provided!'))
        indirCollapsing =ifelse( length(grep('indirCollapsing', args) )>0, strsplit(args[grep('indirCollapsing', args)], '=')[[1]][2], stop('No indirCollapsing argument provided!'))
        inAlignmentPattern =ifelse( length(grep('inAlignmentPattern', args) )>0, strsplit(args[grep('inAlignmentPattern', args)], '=')[[1]][2], stop('No inAlignmentPattern argument provided!'))
        inAmplificationPattern =ifelse( length(grep('inAmplificationPattern', args) )>0, strsplit(args[grep('inAmplificationPattern', args)], '=')[[1]][2], stop('No inAmplificationPattern argument provided!'))
        outDir =ifelse( length(grep('outDir', args) )>0, strsplit(args[grep('outDir', args)], '=')[[1]][2], stop('No outDir argument provided!'))
        sortingOut =ifelse( length(grep('sortingOut', args) )>0, strsplit(args[grep('sortingOut', args)], '=')[[1]][2], stop('No sortingOut argument provided!'))
        amplificationOut =ifelse( length(grep('amplificationOut', args) )>0, strsplit(args[grep('amplificationOut', args)], '=')[[1]][2], stop('No amplificationOut argument provided!'))
        summaryOut =ifelse( length(grep('summaryOut', args) )>0, strsplit(args[grep('summaryOut', args)], '=')[[1]][2], stop('No summaryOut argument provided!'))
}


####
#### ALIGNMENT & OTHER INFO
####

myFiles = list.files(path= indirAlignment, pattern= inAlignmentPattern)

t=as.data.frame(matrix(nrow=length(myFiles), ncol=8))
rownames(t)=do.call(cbind, strsplit(myFiles, '_'))[1,]
colnames(t)=c('input', 'aligned', 'unique       ', 'unique (%)', 'unique human', 'human (% of unique)', 'unique IVT', 'IVT (% of unique)')
for( f in myFiles) {
        lines = readLines(con=paste(indirAlignment,f,sep=''))
        vals=as.numeric(unlist(lapply(strsplit(lines[c(1:3,5:6)], ' '), function(x) x[1])))
        t[strsplit(f, '_')[[1]][1],] = c( vals[1], vals[2], vals[3], round(100*vals[3]/vals[1],digits=2), vals[4], round(100*vals[4]/vals[3],digits=2), vals[5], round(100*vals[5]/vals[3],digits=2) )
}


####
#### PCR AMPLIFICATION AND MOLECULES
####

myFiles = list.files(path= indirCollapsing, pattern= inAmplificationPattern)

tAmp=as.data.frame(matrix(nrow=length(myFiles), ncol=2))
rownames(tAmp)=do.call(cbind, strsplit(myFiles, '_'))[1,]
colnames(tAmp)=c('molecules','positions')


for( f in myFiles) {
        myT = read.delim(paste(indirCollapsing,f,sep=''))
        #molecules=length(myT[,1])
        #positions=length(unique(myT[,1:3])[,1])
        tAmp[strsplit(f, '_')[[1]][1],]=c(length(myT[,1]),length(unique(myT[,1:3])[,1]))
}

t2=cbind(t[,c(1,3)], tAmp)

write.table(t, file=paste(indirAlignment, sortingOut, sep=''), quote=F, sep='\t')
write.table(tAmp, file=paste(outDir, amplificationOut,sep=''),quote=F,sep='\t')
write.table(t2, file=paste(outDir, summaryOut,sep=''),quote=F,sep='\t')



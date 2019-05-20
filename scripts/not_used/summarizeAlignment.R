# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Author: jaerveli
# Date: 21 March 2013
# Description:
#       A simple R script to collect output info from script 'filterAlignment_*' to one file.
#       The output file will be written to dir.
#
# cat /g/steinmetz/project/mTEC_Seq/src/summarizeAlignment.R | R --slave --args
#       dir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/'
#       pattern='uniquelyAligned.info'
#       output='alignmentSummary.info'
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#### PROCESS THE INPUT ARGUMENTS

args = commandArgs(TRUE);

if(length(args)==0){
        stop("No arguments provided!\n")
} else {
        dir =ifelse( length(grep('dir', args) )>0, strsplit(args[grep('dir', args)], '=')[[1]][2], stop('No dir argument provided!'))
        pattern =ifelse( length(grep('pattern', args) )>0, strsplit(args[grep('pattern', args)], '=')[[1]][2], stop('No pattern argument provided!'))
        output =ifelse( length(grep('output', args) )>0, strsplit(args[grep('output', args)], '=')[[1]][2], stop('No output argument provided!'))
}


#### COLLECT THE INFO

myFiles = list.files(path=dir, pattern=pattern)

t=as.data.frame(matrix(nrow=length(myFiles), ncol=8))
rownames(t)=do.call(cbind, strsplit(myFiles, '_'))[1,]
colnames(t)=c('input', 'aligned', 'unique       ', 'unique (%)', 'unique genome', 'genome (% of unique)', 'unique IVT', 'IVT (% of unique)')
for( f in myFiles) {
        lines = readLines(con=paste(dir,f,sep=''))
        vals=as.numeric(unlist(lapply(strsplit(lines[c(1:3,5:6)], ' '), function(x) x[1])))
        t[strsplit(f, '_')[[1]][1],] = c( vals[1], vals[2], vals[3], round(100*vals[3]/vals[1],digits=2), vals[4], round(100*vals[4]/vals[3],digits=2), vals[5], round(100*vals[5]/vals[3],digits=2) )
}

# Write out.
write.table( t, file=paste(dir, output, sep=''), quote=F, sep='\t')

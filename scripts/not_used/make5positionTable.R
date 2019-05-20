# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Author: jaerveli
# Date: 23 March 2013
# Description:
#      A simple script to merge collapsed tables (countsPerMol.tsv) of individual samples
#	to one common table.
#
# cat /g/steinmetz/project/mTEC_Seq/src/make5positionTable.R | R --slave --args
#       indir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/'
#       pattern='_collapsed_countsPerMol.tsv'
#       outdir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/'
#       outfile='collapsed-counts-table_genomic-chr_single-nt.tsv'
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

options(scipen=20)


####### PROCESS THE INPUT ARGUMENTS

args = commandArgs(TRUE);

if(length(args)==0){
        stop("No arguments provided!\n")
} else {
        indir =ifelse( length(grep('indir', args) )>0, strsplit(args[grep('indir', args)], '=')[[1]][2], stop('No indir argument provided!'))
        pattern =ifelse( length(grep('pattern', args) )>0, strsplit(args[grep('pattern', args)], '=')[[1]][2], stop('No pattern argument provided!'))
        outdir =ifelse( length(grep('outdir', args) )>0, strsplit(args[grep('outdir', args)], '=')[[1]][2], stop('No outdir argument provided!'))
        outfile =ifelse( length(grep('outfile', args) )>0, strsplit(args[grep('outfile', args)], '=')[[1]][2], stop('No outfile argument provided!'))

}

####### FIXED PARAMS

genomic=c( '1','10','11','12','13','14','15','16','17','18','19','2','3','4','5','6','7','8','9','20','21','22','X','Y' )


####### READ THE INPUT TABLES AND MAKE ONE BIG OUTPUT TABLE

myFiles = list.files(path=indir, pattern=pattern)

listOfTables = lapply(myFiles, function(f) {
        myT=read.delim(paste(indir,f,sep=''))
        myT = myT[which(as.character(myT$chr) %in% genomic),]
        pos2=apply( cbind( as.character(myT$chr), as.character(myT$strand), as.character(myT$pos) ), 1, function(x) paste(x, collapse='_')) # apply(myT[,1:3], 1, function(x) paste(x, collapse='_'))
        myT2 = cbind(myT,pos2)  
        countsPerPos = table(myT2$pos2)
        df=as.data.frame(cbind(names(countsPerPos), as.integer(countsPerPos)))
        #names(df)=strsplit(f, pattern)[[1]]
        return(df)
})

sampleNames=unlist( lapply(myFiles, function(f) { strsplit(f, pattern)[[1]] }) )

df= listOfTables[[1]]
colnames(df)[1:2]=c('position', sampleNames[1])

for( i in 2:length(listOfTables) ) {    
        cat( paste('Looking at table ', i,'\n', sep='') )
        df = merge( df, listOfTables[[i]], by.x=1, by.y=1, all=T)
        colnames(df)[1+i]=sampleNames[i]
}

# colnames(df)[2:length(listOfTables)] = sampleNames


df2=as.data.frame(matrix(nrow=length(df[,1]), ncol=0))
df2$chr=unlist( lapply(as.character(df[,1]), function(x) strsplit(x, '_')[[1]][1] ) )
df2$strand=unlist( lapply(as.character(df[,1]), function(x) strsplit(x, '_')[[1]][2] ) )
df2$pos=unlist( lapply(as.character(df[,1]), function(x) strsplit(x, '_')[[1]][3] ) )

cols = apply( df[,2:length(df[1,])], 2, function(x) { as.integer(as.character(x)) })

df3 = cbind(df2, cols)
df3[is.na(df3)]=0

## Save the table
colnames(df3)=sub('-', '_', colnames(df3))
write.table(df3, file=paste(outdir, outfile, sep=''), quote=F, sep='\t', row.names=F) 



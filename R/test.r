

args = commandArgs(TRUE);                                                                                                                                                                                                                                                   
print(args)

if(length(args)==0){
    stop("No arguments provided!\n")
    } else {
        dir =ifelse( length(grep('dir', args) )>0, strsplit(args[grep('dir', args)], '=')[[1]][2], stop('No dir argument provided!'))
        pattern =ifelse( length(grep('pattern', args) )>0, strsplit(args[grep('pattern', args)], '=')[[1]][2], stop('No pattern argument provided!'))
        output =ifelse( length(grep('output', args) )>0, strsplit(args[grep('output', args)], '=')[[1]][2], stop('No output argument provided!'))
}

'''
Created on 21 March, 2014

@author: jaerveli

Description:
   
   A script that processes a GSNAP alignment and kicks out reads with ambiguous alignments.
   (this is the usual script for filtering out ambiguous reads for non-TIF-seq samples; here just for human
   (/g/steinmetz/genome/Homo_sapiens/37.68/fasta/DNA/))
   
   Note: This version kicks out ALL reads with ambiguous alignment.
   
Usage:

   python /g/steinmetz/project/mTEC_Seq/src/filterAlignment_human-37.68.py
   --samfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-lo.sam
   --outfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-lo_uniquelyAligned.sam
   --outinfo /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-lo_uniquelyAligned.info

 
'''


import sys, re, HTSeq, optparse

###############################################################
#### READ INPUT
###############################################################

parser = optparse.OptionParser()
parser.add_option('-s', '--samfile', action="store", dest="samfile", type="string", nargs = 1)
parser.add_option('-o', '--outfile', action="store", dest="outfile", type="string", nargs = 1)
parser.add_option('-i', '--outinfo', action="store", dest="outinfo", type="string", nargs = 1)
options, args = parser.parse_args()


mouse_chroms = [ '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'MT', 'X', 'Y' ]
ivt_chroms = [ 'pGIBS-LYS', 'pGIBS-PHE', 'pGIBS-THR', 'gi|49175990_c66000-62000', 'gi|49175990_3526000-3530000' ]


###############################################################
#### PROCESS SAM FILE
###############################################################

#ofile = open(outfile, 'w')
ofile = open(options.outfile, 'w')


########################
# HEADER: read and print
i = 0

# for line in open( samfile ):   
for line in open( options.samfile ):
   line = re.sub("\n", "", line)
   if re.match("@", line) != None:
      ofile.write(line + "\n")
   else:
      break


#######################
# READS: 
# Keep only reads with
# unique alignment

readCounter = 0
unambiguous = 0
ambiguous = 0
aligned = 0
unaligned = 0
mouse = 0
ivts = 0

#sam = HTSeq.SAM_Reader( samfile )
sam = HTSeq.SAM_Reader( options.samfile )
b = HTSeq.bundle_multiple_alignments( sam )


#counter = 0
for r in b:
   
   readCounter += 1
   if readCounter % 100000 == 0:
      sys.stdout.write( str( readCounter ) + " reads processed\n" )
   
   # Read maps uniquely   
   if len(r) == 1:
      
      if r[0].aligned:
         
         if mouse_chroms.count( r[0].iv.chrom ) > 0:
            mouse += 1
         #elif pombe_chroms.count( r[0].iv.chrom ) > 0:
         #   pombe += 1
         elif ivt_chroms.count( r[0].iv.chrom ) > 0:
            ivts += 1
         else:
            std.error.write( 'Unrecognized chromosome ' + r[0].iv.chrom  + '\n' )
            sys.exit()
         
         ofile.write( r[0].original_sam_line )
         unambiguous += 1
      else:
         unaligned += 1
   
   # Ambiguous -- read maps to multiple places
   elif len(r) > 1:
      ambiguous += 1
      #for i in range(len(r)):
      #   sys.stdout.write( r[i].original_sam_line )


ofile.close()



#######################################################################################
# Write info

#info = open(outinfo, 'w')
info = open(options.outinfo, 'w')

info.write( str(readCounter) + " total reads\n" )
info.write( str( unambiguous + ambiguous ) + " total aligned reads\n" )
info.write( str( unambiguous ) + " unambiguous reads\n\n" )   
info.write( str( mouse ) + " unambiguous mouse reads\n" )   
info.write( str( ivts ) + " unambiguous IVT reads\n\n" )   
info.close()


# Std out
sys.stdout.write( str(readCounter) + " total reads\n" )
sys.stdout.write( str( unambiguous + ambiguous ) + " total aligned reads\n" )
sys.stdout.write( str( unambiguous ) + " unambiguous reads\n\n" )   
sys.stdout.write( str( mouse ) + " unambiguous mouse reads\n" )   
sys.stdout.write( str( ivts ) + " unambiguous IVT reads\n\n" )   



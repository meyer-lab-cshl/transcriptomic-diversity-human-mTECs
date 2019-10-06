'''
Created on 22 March 2011

@author: jaerveli

Description:
   
   Splits reads in a sam file based on their strand. Only for aligned reads.
   
   PIPELINE SECTION: VISUALIZATION (IGV)
   
Example:

   python /g/steinmetz/project/Tag-seq/src/divideByStrand.py 
      /g/steinmetz/project/Tag-seq/run/62YA0AAXX_TagSeq_mut1_11s002025-1-1_Pelechano_Garcia_s_1/alignment/62YA0AAXX_TagSeq_mut1_11s002025-1-1_AAGCTA_trimmed_allAligned.sam
      /g/steinmetz/project/Tag-seq/run/62YA0AAXX_TagSeq_mut1_11s002025-1-1_Pelechano_Garcia_s_1/alignment/62YA0AAXX_TagSeq_mut1_11s002025-1-1_AAGCTA_trimmed_allAligned_plus.sam
      /g/steinmetz/project/Tag-seq/run/62YA0AAXX_TagSeq_mut1_11s002025-1-1_Pelechano_Garcia_s_1/alignment/62YA0AAXX_TagSeq_mut1_11s002025-1-1_AAGCTA_trimmed_allAligned_minus.sam


sam_filename = "/g/steinmetz/project/Tag-seq/run/62YA0AAXX_TagSeq_mut1_11s002025-1-1_Pelechano_Garcia_s_1/alignment/62YA0AAXX_TagSeq_mut1_11s002025-1-1_AAGCTA_trimmed_allAligned.sam"
ofilePlus = "/g/steinmetz/project/Tag-seq/run/62YA0AAXX_TagSeq_mut1_11s002025-1-1_Pelechano_Garcia_s_1/alignment/62YA0AAXX_TagSeq_mut1_11s002025-1-1_AAGCTA_trimmed_allAligned_plus.sam"
ofileMinus = "/g/steinmetz/project/Tag-seq/run/62YA0AAXX_TagSeq_mut1_11s002025-1-1_Pelechano_Garcia_s_1/alignment/62YA0AAXX_TagSeq_mut1_11s002025-1-1_AAGCTA_trimmed_allAligned_minus.sam"


Updates:

   3 May 2011: There was a bug: minus strand reads were written to plus files
   29 June 2012: Updated the SAM_Alignment line to be compatible with new version of HTSeq (0.5.3p9).

'''


import sys, re
import HTSeq


sam_filename = sys.argv[1]
ofilePlus = sys.argv[2]
ofileMinus = sys.argv[3]


oPlus = open(ofilePlus, "w")
oMinus = open(ofileMinus, "w")

i = 0
for line in open( sam_filename ):
   line = re.sub("\n", "", line)
   if re.match("@", line) != None:
      oPlus.write(line + "\n")
      oMinus.write(line + "\n")
   else:
      #read = HTSeq.SAM_Alignment( line )
      read = HTSeq.SAM_Alignment.from_SAM_line( line )
      if read.aligned:
         if read.iv.strand == "-":
            oMinus.write(line + "\n")
         else:
            oPlus.write(line + "\n")
   i += 1
   if i % 100000 == 0:
      sys.stderr.write( "%d lines processed.\n" % i )


oPlus.close()
oMinus.close()







'''
Created on 22 March 2011

@author: jaerveli, Hannah Meyer

Description:

   Splits reads in a sam file based on their strand. Only for aligned reads.

Example:

   python divideByStrand.py --samfile [] --outprefix []

sam_filename = "/g/steinmetz/project/Tag-seq/run/62YA0AAXX_TagSeq_mut1_11s002025-1-1_Pelechano_Garcia_s_1/alignment/62YA0AAXX_TagSeq_mut1_11s002025-1-1_AAGCTA_trimmed_allAligned.sam"
outprefix = "/g/steinmetz/project/Tag-seq/run/62YA0AAXX_TagSeq_mut1_11s002025-1-1_Pelechano_Garcia_s_1/alignment/62YA0AAXX_TagSeq_mut1_11s002025-1-1_AAGCTA_trimmed_allAligned"


'''
import sys, HTSeq, optparse, re
import pandas as pd

parser = optparse.OptionParser()
parser.add_option('--outprefix', dest = 'outprefix', type="string", nargs = 1)
parser.add_option('--samfile', dest = 'samfile', type="string", nargs = 1)

options, args = parser.parse_args()

required="outprefix samfile".split()
for r in required:
    if options.__dict__[r] is None:
            parser.error("Parameter %s required but not provided"%r)


sam_filename = options.samfile
ofilePlus = "{}.plus.sam".format(options.outprefix)
ofileMinus = "{}.minus.sam".format(options.outprefix)


oPlus = open(ofilePlus, "w")
oMinus = open(ofileMinus, "w")

i = 0
for line in open( sam_filename ):
   line = re.sub("\n", "", line)
   if re.match("@", line) != None:
      oPlus.write(line + "\n")
      oMinus.write(line + "\n")
   else:
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







'''
Created on 2 April, 2012

@author: jaerveli, Hannah Meyer

Description:
      
      A script to fetch information from an info file (created by get-boundary-reads.py)
      for each read based on read name in input sam file, and attach this information
      as an optional field flag.
      
      Needed at least to attach molecular barcode information to sam files after
      alignment of fastq file.
      
Example:
   
      python /g/steinmetz/project/TSES/src/addFeature.py
      --sam /g/steinmetz/project/TSES/human/run/D194KACXX_PET_human_1_12s006722/alignment/liverI_AA.sam
      --info /g/steinmetz/project/TSES/human/run/D194KACXX_PET_human_1_12s006722/trimming/liverI_AA.info
      --infoname pairing_type
      --flagname T
      --out /g/steinmetz/project/TSES/human/run/D194KACXX_PET_human_1_12s006722/alignment/liverI_AA_pairing.sam

      sam = '/g/steinmetz/project/TSES/human/run/D194KACXX_PET_human_1_12s006722/alignment/liverI_AA.sam'
      info = '/g/steinmetz/project/TSES/human/run/D194KACXX_PET_human_1_12s006722/trimming/liverI_AA.info'
      infoname = 'pairing_type'
      flagname = 'T'
      out = '/g/steinmetz/project/TSES/human/run/D194KACXX_PET_human_1_12s006722/alignment/liverI_AA_pairing.sam'

   Recommendation for flag names for consistency:
      M=molecular_barcode
      P=chim_prox
      C=common_seq --> change this to 1?
      D=chim_dist
      S=common_seq_2 --> change this to 2?
      T=pairing_type
      3=p3_direction
      5=p5_direction
      S=bases_softclipped --> if com seq 2 changed
      C=confidence --> if com seq 1 changed

sam = 'alignment/liverI_3p_AA_alignment_paired.sam'
info = 'trimming/liverI_AA.info'
infoname = 'pairing_type'
flagname = 'M'
out = 'alignment/liverI_3p_AA_alignment_paired_molbc.sam'


sam = 'alignment/liverI_3p_AA_paired_molbc_pairingtype.sam'
info = '/g/steinmetz/project/TSES/human/run/D194KACXX_PET_human_1_12s006722/trimming/liverI_AA.info'
infoname = 'confidence'
flagname = 'C'
out = 'alignment/liverI_3p_AA_paired_molbc_pairingtype_confidence.sam'


'''

# import sys, os, optparse, re, itertools, numpy, HTSeq
import sys, os, optparse, re, HTSeq #, numpy

parser = optparse.OptionParser()
parser.add_option('-s', '--sam', action="store", dest="sam", type="string",
        nargs = 1)
parser.add_option('-i', '--info', action="store", dest="info", type="string",
        nargs = 1)
parser.add_option('-n', '--infoname', action="store", dest="infoname",
        choices = ['molecular_bc', 'chim_prox', 'common_seq', 'chim_dist',
            'common_seq_2', 'pairing_type', 'p3_direction', 'p5_direction',
            'confidence' ])
parser.add_option('-f', '--flagname', action="store", dest="flagname",
        type="string", nargs = 1)
parser.add_option('-o', '--out', action="store", dest="out", type="string",
        nargs = 1)


####################
# READ INPUT
####################

options, args = parser.parse_args()


####################
# DEFINE FUNCTIONS
####################

def parsefile(filename):
    file = open(filename,'r')
    line = file.readline()
    hdrs = re.split('\t',re.sub('[\r\n]','',line))
    rv = {}
    for hd in hdrs:
        rv[hd] = []
    line = file.readline()
    while line!='':
        items = re.split('\t',re.sub('[\r\n]','',line))
        for i in range(len(items)):
            rv[hdrs[i]].append(items[i])
        line = file.readline()
    file.close()
    return rv


###########################
## READ THE INFO FILE   ###
###########################

sys.stdout.write('Reading info file...\n')
table = parsefile(options.info)

# build a dictionary of read name and feature of interest.
infodict = {}
for i in range(len(table[options.infoname])):
   if table['read_name'][i] not in infodict:
      infodict[ table['read_name'][i] ] = table[options.infoname][i]
   else:
      sys.stderr.write( 'Read name ' + table['read_name'][i] + ' present multiple times in info file. Exiting...\n' )
      sys.exit()


###########################
## PROCESS THE SAM FILE ###
###########################

sys.stdout.write( 'Processing sam file...\n' )

## Print out header to output file
ofile = open(options.out, 'w')
i = 0
for line in open( options.sam ):
   line = re.sub("\n", "", line)
   if re.match("@", line) != None:
      ofile.write(line + "\n")
   else:
      break

## Add new optional field to samfile
counter = 0
for read in HTSeq.SAM_Reader( options.sam ):

   counter += 1
   if counter % 100000 == 0:
      sys.stdout.write( str(counter) + " reads processed\n" )

   feature = infodict[ read.read.name ]
   original = re.sub("\n", "", read.original_sam_line)

   # Added this 7 Oct 2012: Possible that this feature has no value.
   # In such a case, print out 'None' -- in this way the sam file syntax stays correct
   if len(feature) == 0:
      newLine = original + '\tX' +  options.flagname + ':' + 'Z:' + 'None' + '\n'
   else:
      newLine = original + '\tX' +  options.flagname + ':' + 'Z:' + feature + '\n'
   ofile.write( newLine )

ofile.close()

sys.stdout.write( 'DONE!\n' )





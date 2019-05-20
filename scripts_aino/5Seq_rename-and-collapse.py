'''
Created on 31 Jan, 2014

@author: jaerveli

Description:
   
   This script can be used to
		a) rename reads based on their 5' position and barcode
		b) build a collapsed sam file (collapsing based on 5' position and barcode)
		c) write a table with collapsed positions, barcode, and molecule count for each
   
   Note: Input arguments --out-renamed, --out-collapsed, and --out-table are optional (in principle).
		 I have not tested giving different inputs though, because usually in the pipeline you'd want all of them.
   
Usage:

   python /g/steinmetz/project/TSES/src/5Seq_rename-and-collapse.py
	--samfile /g/steinmetz/project/TSES/mouse/run/5Seq_2013-10-08-C2GY2ACXX/alignment/46C-2_uniquelyAligned_molbc.sam
	--out-renamed /g/steinmetz/project/TSES/mouse/run/5Seq_2013-10-08-C2GY2ACXX/alignment/46C-2_uniquelyAligned_molbc_renamed.sam
	--out-collapsed /g/steinmetz/project/TSES/mouse/run/5Seq_2013-10-08-C2GY2ACXX/collapsing/46C-2_uniquelyAligned_molbc_collapsed.sam
	--out-molecule-count /g/steinmetz/project/TSES/mouse/run/5Seq_2013-10-08-C2GY2ACXX/collapsing/46C-2_collapsed.tsv
	--out-position-count /g/steinmetz/project/TSES/mouse/run/5Seq_2013-10-08-C2GY2ACXX/collapsing/46C-2_collapsed.tsv
	--out-info /g/steinmetz/project/TSES/mouse/run/5Seq_2013-10-08-C2GY2ACXX/collapsing/46C-2_collapsed.info

samfile = '/g/steinmetz/project/TSES/mouse/run/5Seq_2013-10-08-C2GY2ACXX/alignment/46C-2_uniquelyAligned_molbc.sam'
out_renamed = '/g/steinmetz/project/TSES/mouse/run/5Seq_2013-10-08-C2GY2ACXX/alignment/46C-2_uniquelyAligned_molbc_renamed.sam'
out_collapsed = '/g/steinmetz/project/TSES/mouse/run/5Seq_2013-10-08-C2GY2ACXX/collapsing/46C-2_uniquelyAligned_molbc_collapsed.sam'
out_molecule_count = '/g/steinmetz/project/TSES/mouse/run/5Seq_2013-10-08-C2GY2ACXX/collapsing/46C-2_collapsed_countPerMol.tsv'
out_position_count = '/g/steinmetz/project/TSES/mouse/run/5Seq_2013-10-08-C2GY2ACXX/collapsing/46C-2_collapsed_countPerPos.tsv'
out_info = '/g/steinmetz/project/TSES/mouse/run/5Seq_2013-10-08-C2GY2ACXX/collapsing/46C-2_collapsed.info'

'''


import sys, re, HTSeq, optparse

###############################################################
#### READ INPUT
###############################################################

parser = optparse.OptionParser()
parser.add_option('--samfile', action="store", dest="samfile", type="string", nargs = 1)
parser.add_option('--out-renamed', action="store", dest="out_renamed", type="string", nargs = 1, default = None)
parser.add_option('--out-collapsed', action="store", dest="out_collapsed", type="string", nargs = 1, default = None)
parser.add_option('--out-molecule-count', action="store", dest="out_molecule_count", type="string", nargs = 1, default = None)
parser.add_option('--out-position-count', action="store", dest="out_position_count", type="string", nargs = 1, default = None)
parser.add_option('--out-info', action="store", dest="out_info", type="string", nargs = 1)
options, args = parser.parse_args()



###############################################################
#### PROCESS SAM FILE
###############################################################

if options.out_renamed is not None:
	outAll = open(options.out_renamed, 'w')

if options.out_collapsed is not None:
	outCollapsed = open(options.out_collapsed, 'w')


########################
# HEADER: read and print
i = 0

if (options.out_collapsed is not None) or (options.out_renamed is not None):
	for line in open( options.samfile ):
		line = re.sub("\n", "", line)
		if re.match("@", line) != None:
			if options.out_renamed is not None:
				outAll.write(line + "\n")
			if options.out_collapsed is not None:
				outCollapsed.write(line + "\n")
		else:
			break


#######################
# PROCESS FILE: 
# 

readCounter = 0
nUniqueReads = 0
nUniquePositions = 0

uniqueReads = {}
uniquePositions = {}

sam = HTSeq.SAM_Reader( options.samfile )

#counter = 0
for read in sam:
	#
	readCounter += 1
	if readCounter % 100000 == 0:
		sys.stdout.write( str( readCounter ) + " reads processed\n" )
	#
	## Make new name:
	#  If an gi IVT, clean up the name a little:
	if read.iv.chrom in ['gi|49175990_3526000-3530000', 'gi|49175990_c66000-62000']:
		read.iv.chrom = re.sub( '_', '-', read.iv.chrom )
		read.iv.chrom = re.sub( '\|', '-', read.iv.chrom )
	#
	# Get collapse position
	molbc = read.optional_field('XM')
	if read.iv.strand == '+':
		# Note: + 1 was added here afterwards.
		newName = '_'.join( ( read.iv.chrom, read.iv.strand, str(read.iv.start + 1), molbc ) ) 
	else:
		newName = '_'.join( ( read.iv.chrom, read.iv.strand, str(read.iv.end), molbc ) )
	#
	newSamLine = read.original_sam_line.split('\t')
	newSamLine[0] = newName
	# 
	if options.out_renamed is not None:
		outAll.write( '\t'.join( newSamLine ) )
	#
	## Collect info on unique reads and positions
	## If making collapsed output, write out first encountered instance of each position and BC
	if not uniqueReads.has_key( newName ):
		nUniqueReads += 1
		uniqueReads[ newName ] = 1
		# Write out
		if options.out_collapsed is not None:
			outCollapsed.write( '\t'.join( newSamLine ) )
	else:
		uniqueReads[ newName ] += 1
	# 
	if not uniquePositions.has_key( '_'.join( newName.split('_')[:-1]) ):
		nUniquePositions += 1
		uniquePositions[ '_'.join( newName.split('_')[:-1]) ] = 1
	else:
		uniquePositions[ '_'.join( newName.split('_')[:-1]) ] += 1



#######################
# OUTPUT TABLE: 
# If writing the table, 
# go through the uniqueReads dict.

if options.out_molecule_count is not None:
	outT = open(options.out_molecule_count, 'w')
	outT.write( 'chr\tstrand\tpos\tmolbc\tcount\n' )
	for newName in sorted( uniqueReads.keys() ):
		outT.write( '\t'.join( ( newName.split('_')[0] , newName.split('_')[1], newName.split('_')[2], newName.split('_')[3], str(uniqueReads[ newName ]) ) ) + '\n' )
	outT.close()



#######################
# POSITION SUM TABLE: 

if options.out_position_count is not None:
	outT = open(options.out_position_count, 'w')
	outT.write( 'chr\tstrand\tpos\tcount\n' )
	for newName in sorted( uniquePositions.keys() ):
		outT.write( '\t'.join( ( newName.split('_')[0] , newName.split('_')[1], newName.split('_')[2], str(uniquePositions[ newName ]) ) ) + '\n' )
	outT.close()




#######################
# OUTPUT INFO
# 

info = open(options.out_info, 'w')
#info = open(options.outinfo, 'w')
info.write( str(readCounter) + " total reads\n" )
info.write( str( nUniqueReads ) + " unique molecules\n" )
info.write( str( nUniquePositions ) + " unique positions\n" )
info.close()


# Std out
sys.stdout.write( str(readCounter) + " total reads\n" )
sys.stdout.write( str( nUniqueReads ) + " unique molecules\n" )
sys.stdout.write( str( nUniquePositions ) + " unique positions\n" )




# 
# Date: 30 Jan 2014
# Author: jaerveli
# adjusted: 31.08.2016, Leonie 
#
# Description:
#	
#   Trimming for 5' seq from our lab.
#   Reads the 6bp barcode from the beginning of the reverse file.
#	First 8 bases on the FWD side specify a molecular barcode.
#   
#   Note: 1 mismatch in barcode allowed (if this does not result in ambiguous sample identification).
#	Note: Sample barcode and molecular barcode should be on different sides (fwd and rev reads).
#
#   Input:
#       --forward: raw forward reads from the sequencer
#       --reverse: raw reverse reads from the sequencer
#       --design: Desing file with columns: barcode outfile
#		--sample-barcode: length of sample barcode
#		--molecular-barcode: length of molecular barcode
#		--trim-molecular-barcode: set this if molecular barcode should be trimmed off the read.
#		--info: File were run info will be printed (incl. how many reads per sample.)
#
#   The wetlab protocol is essentially this:
#        P5-XXXXXX-NNNNNNNNNN...NNNNNNNN-A-YYYYYY-P7
#           - Y is the sample barcode
#           - X is the molecular barcode
#           - sequencing is paired end.
#	
# Example
#
#	python /g/steinmetz/project/mTEC_Seq/src/5seq_sort-and-trim.py
#       --forward /g/steinmetz/incoming/solexa/2014-03-18-C3V0VACXX/C3V0VACXX_5TSeq_mTec1_13s007938-1-1_Clauder-Muenster_lane313s007938_1_sequence.txt.gz
#       --reverse /g/steinmetz/incoming/solexa/2014-03-18-C3V0VACXX/C3V0VACXX_5TSeq_mTec1_13s007938-1-1_Clauder-Muenster_lane313s007938_2_sequence.txt.gz
#       --design /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/mTEC_5Seq_2014-03-18.design
#		--sample-barcode 6
#		--molecular-barcode 8
#		--trim-molecular-barcode
#		--info /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/mTEC_5Seq_2014-03-18.info
#
###########################

import sys, HTSeq, optparse, re, itertools, difflib

parser = optparse.OptionParser()
parser.add_option('--forward', dest = 'forward', type="string", nargs = 1)
parser.add_option('--reverse', dest = 'reverse', type="string", nargs = 1)
parser.add_option('--design', dest = 'design', type="string", nargs = 1)
parser.add_option('--sample-barcode', dest = 'sample_barcode', type="int", nargs = 1, default = 6)
parser.add_option('--molecular-barcode', dest = 'molecular_barcode', type="int", nargs = 1, default = 8)
parser.add_option('--trim-molecular-barcode', dest='trim_molecular_barcode', action="store_true", default='False')
parser.add_option('--info', dest = 'info', type="string", nargs = 1)

options, args = parser.parse_args()

## DEFINE FUNCTIONS

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

table = parsefile(options.design)
#table = parsefile(design)

forward_file = HTSeq.FastqReader( options.forward )
reverse_file = HTSeq.FastqReader( options.reverse )
#forward_file = HTSeq.FastqReader( forward )
#reverse_file = HTSeq.FastqReader( reverse )

#dictionaries for forward and reverse read files
ofiles_fwd = {}
ofiles_rev = {}
bcfiles = {}
for i in range(len(table['barcode'])):
   ofiles_fwd[table['barcode'][i]] = open( table['outfile_fwd'][i], "w" )
   ofiles_rev[table['barcode'][i]] = open( table['outfile_rev'][i], "w" )
   bcfiles[table['barcode'][i]] = open( table['molecular_bc_file'][i], "w" )

# Write header to molecular barcode files
# Note: Added 22 March 2014 AFTER running the script
# in order to correct the output to molecular barcode file.
# NOT TESTED!
#for bcf in bcfiles:

#added by Leonie: iterate over values in dict instead of over keys
for bcf in bcfiles.values():
   bcf.write('read_name\tmolecular_bc\n')
   #ofiles[table['barcode'][i]] = open( table['outfile'][i], "w" )
   #bcfiles[table['barcod'][i]] = open( table['molecular_bc_file'][i], "w" )

# Go through the reads and split to output files according to barcode

counter = 0
exact_match = 0
one_mismatch = 0
no_match = 0
readsPerSample = {}

for fwd, rev in itertools.izip( forward_file, reverse_file ):
	counter += 1
	if counter % 100000 == 0:
		sys.stderr.write(str(counter) + ' reads processed\n')
		#break
	#
	## Confirm that reads are pairs
	# 
	rfw_name = fwd.name.split(" ")[0]
	rrev_name = rev.name.split(" ")[0]
	assert rfw_name == rrev_name
	#
	## EXACT MATCH barcode on REV. Write out FWD.
	#
	if ofiles_fwd.has_key( rev.seq[0:options.sample_barcode] ):
		exact_match += 1
		if options.trim_molecular_barcode:
			fwd[options.molecular_barcode:].write_to_fastq_file( ofiles_fwd[ rev.seq[0:options.sample_barcode] ] )
			rev[7:].write_to_fastq_file( ofiles_rev[ rev.seq[0:options.sample_barcode] ] )
		else:
			fwd.write_to_fastq_file( ofiles_fwd[ rev.seq[0:options.sample_barcode] ] )
			rev[7:].write_to_fastq_file( ofiles_rev[ rev.seq[0:options.sample_barcode] ] )
		# 
		## Write out information about the molecular barcode 
		# 
		#bcfiles[ rev.seq[0:options.sample_barcode] ].write( fwd.name + '\t' + fwd.seq[0:options.molecular_barcode] + '\n' )
		# Note: Added 22 March 2014 AFTER running the script in order to correct the output to molecular barcode file.
		bcfiles[ rev.seq[0:options.sample_barcode] ].write( rfw_name + '\t' + fwd.seq[0:options.molecular_barcode] + '\n' )
		#
		## Keep count of reads per sample:
		# 
		if readsPerSample.has_key( rev.seq[0:options.sample_barcode] ):
			readsPerSample[ rev.seq[0:options.sample_barcode] ] += 1
		else:
			readsPerSample[ rev.seq[0:options.sample_barcode] ] = 1
	#
	## NO EXACT MATCH barcode on REV. Write out FWD.
	#
	else:
		matchWithMismatch = 0
		bcMatched = ''
		# 
		for i in range(len(table['barcode'])):
			if difflib.SequenceMatcher(None, table['barcode'][i], rev.seq[0:options.sample_barcode]).ratio() >= 0.80: # Allowed ratio of mismatch.
				matchWithMismatch += 1
				bcMatched=table['barcode'][i]
		# 
		if matchWithMismatch == 1:
			one_mismatch += 1

			if options.trim_molecular_barcode:
				fwd[options.molecular_barcode:].write_to_fastq_file( ofiles_fwd[ bcMatched ] )
				rev[7:].write_to_fastq_file( ofiles_rev[ bcMatched ] )
			else:
				fwd.write_to_fastq_file( ofiles_fwd[ bcMatched ] )
				rev[7:].write_to_fastq_file( ofiles_rev[ bcMatched ] )
			# 
			## Write out information about the molecular barcode 
			# 
			# bcfiles[ rev.seq[0:sample_barcode] ].write( fwd.name + '\t' + rev.seq[0:molecular_barcode] + '\n' )
			#bcfiles[ bcMatched ].write( fwd.name + '\t' + fwd.seq[0:options.molecular_barcode] + '\n' )
			# Note: Added 22 March 2014 AFTER running the script in order to correct the output to molecular barcode file.
			bcfiles[ bcMatched ].write( rfw_name + '\t' + fwd.seq[0:options.molecular_barcode] + '\n' )
			#
			## Keep count of reads per sample:
			# 
			if readsPerSample.has_key( bcMatched ):
				readsPerSample[ bcMatched ] += 1
			else:
				readsPerSample[ bcMatched ] = 1
		# 
		else:
			no_match += 1

### Close output files
for f in ofiles_fwd.values():
   f.close()

for f in ofiles_rev.values():
   f.close()

for f in bcfiles.values():
   f.close()



### Write out info

# Into an info file
oinfo = open(options.info, 'w')
#oinfo = open(info, 'w')
oinfo.write( str(counter) + ' reads in total' + '\n' )
oinfo.write( str(exact_match) + ' reads with exact match to barcode' + '\n' )
oinfo.write( str(one_mismatch) + ' reads with non-exact match to barcode (one mismatch allowed)' + '\n' )
oinfo.write( str(no_match) + ' reads without a match to barcode' + '\n\n' )

oinfo.write( 'barcode\tsample\n' )
for sample in sorted(readsPerSample.items(), key=lambda x: x[1]):
	oinfo.write( sample[0] + '\t' + str(sample[1]) + '\n' )

oinfo.close()



# Into standard out.
sys.stdout.write( str(counter) + ' reads in total' + '\n' )
sys.stdout.write( str(exact_match) + ' reads with exact match to barcode' + '\n' )
sys.stdout.write( str(one_mismatch) + ' reads with non-exact match to barcode (one mismatch allowed)' + '\n' )
sys.stdout.write( str(no_match) + ' reads without a match to barcode' + '\n' )

sys.stdout.write( 'barcode\tsample\n' )
for sample in sorted(readsPerSample.items(), key=lambda x: x[1]):
	sys.stdout.write( sample[0] + '\t' + str(sample[1]) + '\n' )


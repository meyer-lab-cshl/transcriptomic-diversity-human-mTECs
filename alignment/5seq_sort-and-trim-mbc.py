# Author: jaerveli,Leonie Stroemich,Hannah Meyer
#

# Description:
#   Trimming for 5' seq data
#   Reads the 6bp barcode from the beginning of the reverse file.
#   First 8 bases on the FWD side specify a molecular barcode.
#   Note: 1 mismatch in barcode allowed (if this does not result in ambiguous sample identification).
#       Note: Sample barcode and molecular barcode should be on different sides (fwd and rev reads).
#
#   Input:
#       --forward: raw forward reads from the sequencer
#       --reverse: raw reverse reads from the sequencer
#       --multiplexed: boolean, indicating if reads are from multiplexed experiment
#       --design: Tab-sep design file with col: barcode outfile_fwd outfile_rev
#       --sample-barcode: length of sample barcode
#       --molecular-barcode: length of molecular barcode
#       --trim-molecular-barcode: set this if molecular barcode should be trimmed off the read.
#       --info: File were run info will be printed (incl. how many reads per sample.)
#   The wetlab protocol is essentially this:
#        P5-XXXXXX-NNNNNNNNNN...NNNNNNNN-A-YYYYYY-P7
#           - Y is the sample barcode
#           - X is the molecular barcode
#           - sequencing is paired end.
# Example
#
#       python /g/steinmetz/project/mTEC_Seq/src/5seq_sort-and-trim.py
#           --forward /g/steinmetz/incoming/solexa/2014-03-18-C3V0VACXX/C3V0VACXX_5TSeq_mTec1_13s007938-1-1_Clauder-Muenster_lane313s007938_1_sequence.txt.gz
#           --reverse /g/steinmetz/incoming/solexa/2014-03-18-C3V0VACXX/C3V0VACXX_5TSeq_mTec1_13s007938-1-1_Clauder-Muenster_lane313s007938_2_sequence.txt.gz
#           --design /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/mTEC_5Seq_2014-03-18.design
#           --sample-barcode 6
#           --molecular-barcode 8
#           --trim-molecular-barcode
#           --info /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/mTEC_5Seq_2014-03-18.info
#
###########################

import sys, HTSeq, optparse, re
import difflib
import pdb
import pandas as pd

parser = optparse.OptionParser()
parser.add_option('--forward', dest = 'forward', type="string", nargs = 1)
parser.add_option('--reverse', dest = 'reverse', type="string", nargs = 1)
parser.add_option('--design', dest = 'design', type="string", nargs = 1)
parser.add_option('--sample-barcode', dest = 'sample_barcode', type="int",
        nargs = 1, default = 6)
parser.add_option('--molecular-barcode', dest = 'molecular_barcode', type="int",
        nargs = 1, default = 8)
parser.add_option('--trim-molecular-barcode', dest='trim_molecular_barcode',
        action="store_true", default=False)
parser.add_option('--multiplexed', dest='multiplexed',
        action="store_true", default=False)
parser.add_option('--info', dest = 'info', type="string", nargs = 1)

options, args = parser.parse_args()

table = pd.read_csv(options.design, sep=" ")
forward_file = HTSeq.FastqReader(options.forward)
reverse_file = HTSeq.FastqReader(options.reverse)
if options.multiplexed:
    ofiles_fwd = {}
    ofiles_rev = {}
    bcfiles = {}
    for i in range(len(table['barcode'])):
        ofiles_fwd[table['barcode'][i]] = open(table['outfile_fwd'][i], "w")
        ofiles_rev[table['barcode'][i]] = open(table['outfile_rev'][i], "w")
        bcfiles[table['barcode'][i]] = open(table['molecular_bc_file'][i], "w")
    # Write header to molecular barcode files
    for bcf in bcfiles.values():
        bcf.write('read_name\tmolecular_bc\n')
else:
    ofiles_fwd = open(table['outfile_fwd'][0], "w")
    ofiles_rev = open(table['outfile_rev'][0], "w")
    bcfiles = open(table['molecular_bc_file'][0], "w")
    bcfiles.write('read_name\tmolecular_bc\n')


# Go through the reads and split to output files according to barcode
ratio = (options.sample_barcode -1)/options.sample_barcode
counter = 0
exact_match = 0
one_mismatch = 0
no_match = 0
readsPerSample = {}

for fwd, rev in zip(forward_file, reverse_file):
    counter += 1
    if counter % 100000 == 0:
        sys.stderr.write(str(counter) + ' reads processed\n')

    ## Confirm that reads are pairs
    rfw_name = fwd.name.split(" ")[0]
    rrev_name = rev.name.split(" ")[0]
    assert rfw_name == rrev_name

    bc_mol = fwd.seq.decode('utf-8')[0:options.molecular_barcode]
    if options.multiplexed :
        ## beginnings of reads correspond to sample (rev) and molecular (fwd) bc
        bc_sample = rev.seq.decode('utf-8')[0:options.sample_barcode]
        if bc_sample in ofiles_rev.keys():
            exact_match += 1
            rev.write_to_fastq_file(ofiles_rev[bc_sample])
            if options.trim_molecular_barcode:
                fwd[options.molecular_barcode:].write_to_fastq_file(
                        ofiles_fwd[bc_sample])
            else:
                fwd.write_to_fastq_file(ofiles_fwd[bc_sample])
            rev.write_to_fastq_file(ofiles_rev[bc_sample])

            ## Write out information about the molecular barcode
            bcfiles[bc_sample].write(rfw_name + '\t' + bc_mol + '\n' )

            ## Keep count of reads per sample:
            if bc_sample in readsPerSample.keys():
                readsPerSample[bc_sample] += 1
            else:
                readsPerSample[bc_sample] = 1
        else:
            matchWithMismatch = 0
            bc_matched = ''
            for i in range(len(table['barcode'])):
                if difflib.SequenceMatcher(None, table['barcode'][i],
                        bc_sample).ratio() >= ratio:
                    matchWithMismatch += 1
                    bc_matched=table['barcode'][i]

            if matchWithMismatch == 1:
                one_mismatch += 1
                if options.trim_molecular_barcode:
                    fwd[options.molecular_barcode:].write_to_fastq_file(
                            ofiles_fwd[bc_matched])
                else:
                    fwd.write_to_fastq_file(ofiles_fwd[bc_matched])

                ## Write out information about the molecular barcode 
                bcfiles[bc_matched].write(rfw_name + '\t' + bc_mol + '\n' )

                ## Keep count of reads per sample:
                if bc_matched in readsPerSample.keys():
                    readsPerSample[bc_matched] += 1
                else:
                    readsPerSample[bc_matched] = 1
            else:
                no_match += 1
                bc_sample = 'nomatch'
                if options.trim_molecular_barcode:
                    fwd[options.molecular_barcode:].write_to_fastq_file(
                            ofiles_fwd[bc_sample])
                else:
                    fwd.write_to_fastq_file(ofiles_fwd[bc_sample])
                rev.write_to_fastq_file(ofiles_rev[bc_sample])

                ## Write out information about the molecular barcode
                bcfiles[bc_sample].write(rfw_name + '\t' + bc_mol + '\n' )
    else:
        if options.trim_molecular_barcode:
            fwd[options.molecular_barcode:].write_to_fastq_file(
                    ofiles_fwd)
        else:
            fwd.write_to_fastq_file(ofiles_fwd)
        rev.write_to_fastq_file(ofiles_rev)

        ## Write out information about the molecular barcode
        bcfiles.write(rfw_name + '\t' + bc_mol + '\n' )

### Close output files
for f in ofiles_fwd.values():
   f.close()
for f in ofiles_rev.values():
   f.close()
for f in bcfiles.values():
   f.close()

### Write out info

if options.multiplexed:
    # Into an info file
    oinfo = open(options.info, 'w')
    oinfo.write(str(counter) + ' reads in total' + '\n' )
    oinfo.write(str(exact_match) + ' reads with exact match to barcode' + '\n' )
    oinfo.write(str(one_mismatch) +
        ' reads with non-exact match to barcode (one mismatch allowed)' + '\n' )
    oinfo.write(str(no_match) + ' reads without a match to barcode' + '\n\n' )

    oinfo.write('barcode\tsample\n' )
    for sample in sorted(readsPerSample.items(), key=lambda x: x[1]):
        oinfo.write(sample[0] + '\t' + str(sample[1]) + '\n' )
    oinfo.close()

    # Into standard out.
    sys.stdout.write(str(counter) + ' reads in total' + '\n' )
    sys.stdout.write(str(exact_match) +
            ' reads with exact match to barcode' + '\n' )
    sys.stdout.write( str(one_mismatch) +
        ' reads with non-exact match to barcode (one mismatch allowed)' + '\n' )
    sys.stdout.write(str(no_match) + ' reads without a match to barcode' + '\n' )

    sys.stdout.write('barcode\tsample\n' )
    for sample in sorted(readsPerSample.items(), key=lambda x: x[1]):
        sys.stdout.write( sample[0] + '\t' + str(sample[1]) + '\n' )
else:
    # Into an info file
    oinfo = open(options.info, 'w')
    oinfo.write(str(counter) + ' reads in total' + '\n' )
    oinfo.close()
    sys.stdout.write(str(counter) + ' reads in total' + '\n' )


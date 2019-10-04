


############## 
############## GENERAL INFO
############## 

## Description:

High and low AIRE samples of three individuals for Philip's mTEC project.
Sequenced on one HiSeq lane (run read on 2014-March-18), 6-plex.


## Project folder(s):

/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX
/g/tier2/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX


## Raw data:

/g/steinmetz/incoming/solexa/2014-03-18-C3V0VACXX/C3V0VACXX_5TSeq_mTec1_13s007938-1-1_Clauder-Muenster_lane313s007938_1_sequence.txt.gz
/g/steinmetz/incoming/solexa/2014-03-18-C3V0VACXX/C3V0VACXX_5TSeq_mTec1_13s007938-1-1_Clauder-Muenster_lane313s007938_2_sequence.txt.gz


## Sample description on lab wiki:

http://steinmetzlab.embl.de/wiki/index.php/5TSeq_mTec

11	87-lo	11	TTAATT	0.45	5.56
12	87-hi	12	CCTCCC	0.7	3.57
13	212-lo	13	TATATC	0.96	2.6
14	212-hi	14	TGCCGA	1.39	1.8
15	214-lo	15	TGACAT	0.51	4.9
16	214-hi	16	CGCCTG	0.765	3.27






############## 
############## QUALITIES
############## 

cd /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/qualities/

nohup htseq-qa -m 41 -t fastq /g/steinmetz/incoming/solexa/2014-03-18-C3V0VACXX/C3V0VACXX_5TSeq_mTec1_13s007938-1-1_Clauder-Muenster_lane313s007938_1_sequence.txt.gz > nohup.out.fwd &
nohup htseq-qa -m 41 -t fastq /g/steinmetz/incoming/solexa/2014-03-18-C3V0VACXX/C3V0VACXX_5TSeq_mTec1_13s007938-1-1_Clauder-Muenster_lane313s007938_2_sequence.txt.gz > nohup.out.rev &






############## 
############## PRIMED-DIMER
############## 

cd /g/steinmetz/project/TSES/mouse/run/5Seq_2013-10-08-C2GY2ACXX/

zgrep "^@HWI-ST999:171:C3V0VACXX" /g/steinmetz/incoming/solexa/2014-03-18-C3V0VACXX/C3V0VACXX_5TSeq_mTec1_13s007938-1-1_Clauder-Muenster_lane313s007938_1_sequence.txt.gz | wc -l > primerDimer/fwd.readCount.txt

zgrep "GCAGGAATGCCGAG" /g/steinmetz/incoming/solexa/2014-03-18-C3V0VACXX/C3V0VACXX_5TSeq_mTec1_13s007938-1-1_Clauder-Muenster_lane313s007938_1_sequence.txt.gz | wc -l > primerDimer/fwd.primerDimerCount.txt
1914013


100*(1914013/157669770)
[1] 1.213938




############## 
############## SORTING
############## 

## BC at beginning of reverse. Position 4 has huge amount of N. Allow 1 mismatch OR ignore position 1.
## Molecular barcode (8 nt) at the beginning of the FWD. Write that to a separate file: read_name molecular_bc
## All that is already done by an existing script '/g/steinmetz/project/TSES/src/5seq_sort-and-trim.py'
## Copy that to mTEC folder /g/steinmetz/project/mTEC_Seq/src/


## This is the sample design:

cat > /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/mTEC_5Seq_2014-03-18.design

barcode	outfile	molecular_bc_file
TTAATT	/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/87-lo.fastq	/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/87-lo.molbc.tsv
CCTCCC	/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/87-hi.fastq	/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/87-hi.molbc.tsv
TATATC	/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/212-lo.fastq	/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/212-lo.molbc.tsv
TGCCGA	/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/212-hi.fastq	/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/212-hi.molbc.tsv
TGACAT	/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/214-lo.fastq	/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/214-lo.molbc.tsv
CGCCTG	/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/214-hi.fastq	/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/214-hi.molbc.tsv


## Go to project folder:
cd /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX


## Run 5seq_sort-and-trim.py with nohup
### (sorting this big raw data file takes some time (several hours), so you want to be able to log out of the server without interrupting your job)
nohup python /g/steinmetz/project/mTEC_Seq/src/5seq_sort-and-trim.py --forward /g/steinmetz/incoming/solexa/2014-03-18-C3V0VACXX/C3V0VACXX_5TSeq_mTec1_13s007938-1-1_Clauder-Muenster_lane313s007938_1_sequence.txt.gz --reverse /g/steinmetz/incoming/solexa/2014-03-18-C3V0VACXX/C3V0VACXX_5TSeq_mTec1_13s007938-1-1_Clauder-Muenster_lane313s007938_2_sequence.txt.gz --design /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/mTEC_5Seq_2014-03-18.design --sample-barcode 6 --molecular-barcode 8 --trim-molecular-barcode --info /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/mTEC_5Seq_2014-03-18.info > sorting/nohup.mTEC_5Seq_2014-03-18.out &


157669770 reads in total
102315606 reads with exact match to barcode
36882975 reads with non-exact match to barcode (one mismatch allowed)
18471189 reads without a match to barcode

barcode sample
TTAATT  14420646
CGCCTG  21057102
CCTCCC  23116032
TGACAT  23554099
TATATC  23596949
TGCCGA  33453753


cat > /g/steinmetz/project/mTEC_Seq/src/fixMolBCfile.py

import sys, optparse, re

parser = optparse.OptionParser()
parser.add_option('--infile', dest = 'infile', type="string", nargs = 1)
parser.add_option('--outfile', dest = 'outfile', type="string", nargs = 1)

options, args = parser.parse_args()

infile = open(options.infile,'r')
outfile=open(options.outfile,'w')
outfile.write('read_name\tmolecular_bc\n')

line = infile.readline()
while line!='':
	name=re.split('\t', re.sub('[\r\n]','',line))[0]
	name2=re.split(' ', name)[0]
	bc=re.split('\t', re.sub('[\r\n]','',line))[1]
	outfile.write( name2 + '\t' + bc + '\n' )
	line = infile.readline()

outfile.close()
infile.close()


## There was a problem in formatting of the molecular barcode output files.
## That is now fixed in the script.
## This script fixes the output directly so that there is no need to rerun the sortin
python /g/steinmetz/project/mTEC_Seq/src/fixMolBCfile.py --infile sorting/212-hi.molbc.tsv --outfile sorting/212-hi.molbc.2.tsv
python /g/steinmetz/project/mTEC_Seq/src/fixMolBCfile.py --infile sorting/212-lo.molbc.tsv --outfile sorting/212-lo.molbc.2.tsv

python /g/steinmetz/project/mTEC_Seq/src/fixMolBCfile.py --infile sorting/214-hi.molbc.tsv --outfile sorting/214-hi.molbc.2.tsv
python /g/steinmetz/project/mTEC_Seq/src/fixMolBCfile.py --infile sorting/214-lo.molbc.tsv --outfile sorting/214-lo.molbc.2.tsv

python /g/steinmetz/project/mTEC_Seq/src/fixMolBCfile.py --infile sorting/87-hi.molbc.tsv --outfile sorting/87-hi.molbc.2.tsv
python /g/steinmetz/project/mTEC_Seq/src/fixMolBCfile.py --infile sorting/87-lo.molbc.tsv --outfile sorting/87-lo.molbc.2.tsv




############## 
############## QUALITIES FOR SORTED
############## 

cd /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/qualities/

htseq-qa -m 41 -t fastq /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/87-lo.fastq.gz
htseq-qa -m 41 -t fastq /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/87-hi.fastq.gz
htseq-qa -m 41 -t fastq /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/212-lo.fastq.gz
htseq-qa -m 41 -t fastq /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/212-hi.fastq.gz
htseq-qa -m 41 -t fastq /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/214-lo.fastq.gz
htseq-qa -m 41 -t fastq /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/214-hi.fastq.gz


############## 
############## ALIGNMENT
############## 

#######
####### RUN THE ALIGNMENT
#######

# Run on submaster using 50 cores each run.

cat > /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/src/run-alignment_5Seq_mTEC.sh

*~*~*~*

# Author: jaerveli
# Date: 21 March 2014
# Description: Running GSNAP on Philip's mTEC data on submaster with 50 cores per run.

#! /bin/csh -f

cd /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/

### 87-lo
/g/huber/users/jgehring/software/gmap/bin-sch/gmap-2012-01-11/bin/gsnap -t 50 --quality-protocol=sanger -A sam -m 0.07 -B 5 -n 2 -D /g/steinmetz/genome/Homo_sapiens/37.68/indexes/gsnap/Homo_sapiens.GRCh37.68_ivts -d Homo_sapiens.GRCh37.68_ivts -s /g/steinmetz/genome/Homo_sapiens/37.68/indexes/gsnap/Homo_sapiens.GRCh37.68_ivts/Homo_sapiens_37.68_splice-sites.iit sorting/87-lo.fastq > alignment/87-lo.sam

### 87-hi
/g/huber/users/jgehring/software/gmap/bin-sch/gmap-2012-01-11/bin/gsnap -t 50 --quality-protocol=sanger -A sam -m 0.07 -B 5 -n 2 -D /g/steinmetz/genome/Homo_sapiens/37.68/indexes/gsnap/Homo_sapiens.GRCh37.68_ivts -d Homo_sapiens.GRCh37.68_ivts -s /g/steinmetz/genome/Homo_sapiens/37.68/indexes/gsnap/Homo_sapiens.GRCh37.68_ivts/Homo_sapiens_37.68_splice-sites.iit sorting/87-hi.fastq > alignment/87-hi.sam

### 212-lo
/g/huber/users/jgehring/software/gmap/bin-sch/gmap-2012-01-11/bin/gsnap -t 50 --quality-protocol=sanger -A sam -m 0.07 -B 5 -n 2 -D /g/steinmetz/genome/Homo_sapiens/37.68/indexes/gsnap/Homo_sapiens.GRCh37.68_ivts -d Homo_sapiens.GRCh37.68_ivts -s /g/steinmetz/genome/Homo_sapiens/37.68/indexes/gsnap/Homo_sapiens.GRCh37.68_ivts/Homo_sapiens_37.68_splice-sites.iit sorting/212-lo.fastq > alignment/212-lo.sam

### 212-hi
/g/huber/users/jgehring/software/gmap/bin-sch/gmap-2012-01-11/bin/gsnap -t 50 --quality-protocol=sanger -A sam -m 0.07 -B 5 -n 2 -D /g/steinmetz/genome/Homo_sapiens/37.68/indexes/gsnap/Homo_sapiens.GRCh37.68_ivts -d Homo_sapiens.GRCh37.68_ivts -s /g/steinmetz/genome/Homo_sapiens/37.68/indexes/gsnap/Homo_sapiens.GRCh37.68_ivts/Homo_sapiens_37.68_splice-sites.iit sorting/212-hi.fastq > alignment/212-hi.sam

### 214-lo
/g/huber/users/jgehring/software/gmap/bin-sch/gmap-2012-01-11/bin/gsnap -t 50 --quality-protocol=sanger -A sam -m 0.07 -B 5 -n 2 -D /g/steinmetz/genome/Homo_sapiens/37.68/indexes/gsnap/Homo_sapiens.GRCh37.68_ivts -d Homo_sapiens.GRCh37.68_ivts -s /g/steinmetz/genome/Homo_sapiens/37.68/indexes/gsnap/Homo_sapiens.GRCh37.68_ivts/Homo_sapiens_37.68_splice-sites.iit sorting/214-lo.fastq > alignment/214-lo.sam

### 214-hi
/g/huber/users/jgehring/software/gmap/bin-sch/gmap-2012-01-11/bin/gsnap -t 50 --quality-protocol=sanger -A sam -m 0.07 -B 5 -n 2 -D /g/steinmetz/genome/Homo_sapiens/37.68/indexes/gsnap/Homo_sapiens.GRCh37.68_ivts -d Homo_sapiens.GRCh37.68_ivts -s /g/steinmetz/genome/Homo_sapiens/37.68/indexes/gsnap/Homo_sapiens.GRCh37.68_ivts/Homo_sapiens_37.68_splice-sites.iit sorting/214-hi.fastq > alignment/214-hi.sam

*~*~*~*

bsub -M 28000 -R "select[(mem >= 30000)]" -oo /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/src/lsf_output_align-mTEC-5Seq.txt -eo /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/src/lsf_error_align-mTEC-5Seq.txt -n 50 -J "5Seq mTEC" "sh /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/src/run-alignment_5Seq_mTEC.sh"



cat > /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/src/gzip_sorted.molbc.sh
*~*~*~*
#! /bin/csh -f
cd /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/
for i in sorting/*molbc.tsv; do echo $i; gzip $i; done
*~*~*~*
nohup sh /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/src/gzip_sorted.molbc.sh > /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/src/gzip_sorted.molbc.log &
[1] 53799



#######
####### FILTER THE ALIGNMENT (REMOVE AMBIGUOUS READS)
#######


### 87-lo
python /g/steinmetz/project/mTEC_Seq/src/filterAlignment_human-37.68.py --samfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-lo.sam --outfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-lo_uniquelyAligned.sam --outinfo /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-lo_uniquelyAligned.info

14420646 total reads
13768083 total aligned reads
11527255 unambiguous reads

10221313 unambiguous mouse reads
1305942 unambiguous IVT reads

### 87-hi
python /g/steinmetz/project/mTEC_Seq/src/filterAlignment_human-37.68.py --samfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-hi.sam --outfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-hi_uniquelyAligned.sam --outinfo /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-hi_uniquelyAligned.info

23116032 total reads
22262219 total aligned reads
18960575 unambiguous reads

17615214 unambiguous mouse reads
1345361 unambiguous IVT reads

### 212-lo
python /g/steinmetz/project/mTEC_Seq/src/filterAlignment_human-37.68.py --samfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/212-lo.sam --outfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/212-lo_uniquelyAligned.sam --outinfo /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/212-lo_uniquelyAligned.info

23596949 total reads
23011773 total aligned reads
19614602 unambiguous reads

18232011 unambiguous mouse reads
1382591 unambiguous IVT reads

### 212-hi
python /g/steinmetz/project/mTEC_Seq/src/filterAlignment_human-37.68.py --samfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/212-hi.sam --outfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/212-hi_uniquelyAligned.sam --outinfo /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/212-hi_uniquelyAligned.info

33453753 total reads
33011216 total aligned reads
28476242 unambiguous reads

27456942 unambiguous mouse reads
1019300 unambiguous IVT reads


### 214-lo
python /g/steinmetz/project/mTEC_Seq/src/filterAlignment_human-37.68.py --samfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/214-lo.sam --outfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/214-lo_uniquelyAligned.sam --outinfo /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/214-lo_uniquelyAligned.info

23554099 total reads
22414340 total aligned reads
18820555 unambiguous reads

17610074 unambiguous mouse reads
1210481 unambiguous IVT reads

### 214-hi
python /g/steinmetz/project/mTEC_Seq/src/filterAlignment_human-37.68.py --samfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/214-hi.sam --outfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/214-hi_uniquelyAligned.sam --outinfo /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/214-hi_uniquelyAligned.info

21057102 total reads
20169032 total aligned reads
17215709 unambiguous reads

15978501 unambiguous mouse reads
1237208 unambiguous IVT reads


#######
####### GET SUMMARY OF ALIGNMENT
#######

cat /g/steinmetz/project/mTEC_Seq/src/summarizeAlignment.R | R --slave --args dir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/' pattern='uniquelyAligned.info' output='alignmentSummary.info'

212-hi  33453753        33011216        28476242        85.12   27456942        96.42   1019300 3.58
212-lo  23596949        23011773        19614602        83.12   18232011        92.95   1382591 7.05
214-hi  21057102        20169032        17215709        81.76   15978501        92.81   1237208 7.19
214-lo  23554099        22414340        18820555        79.9    17610074        93.57   1210481 6.43
87-hi   23116032        22262219        18960575        82.02   17615214        92.9    1345361 7.1
87-lo   14420646        13768083        11527255        79.94   10221313        88.67   1305942 11.33




############## 
############## DO SOME CLEANING UP
############## 

### GZIP SORTED 

for i in sorting/*.fastq; do echo $i; gzip $i; done


### Store alignment files as bam.

for i in alignment/*hi.sam; do echo $i; samtools view -bS -o alignment/`basename $i .sam`.bam $i; done
for i in alignment/*lo.sam; do echo $i; samtools view -bS -o alignment/`basename $i .sam`.bam $i; done
for i in alignment/*uniquelyAligned.sam; do echo $i; samtools view -bS -o alignment/`basename $i .sam`.bam $i; done

# rm alignment/*hi.sam
# rm alignment/*lo.sam

# Note: Consider removing the non-filtered alignment files (not needed anymore)



############## 
############## ADD MOLECULAR BARCODE
############## 

## Copy the existing script for this purpose from
## /g/steinmetz/project/TSES/src/addFeature.py to python /g/steinmetz/project/mTEC_Seq/src/addFeature.py



### 87-lo
python /g/steinmetz/project/mTEC_Seq/src/addFeature.py --sam /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-lo_uniquelyAligned.sam --info /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/87-lo.molbc.2.tsv --infoname molecular_bc --flagname M --out /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-lo_uniquelyAligned_molbc.sam

### 87-hi
python /g/steinmetz/project/mTEC_Seq/src/addFeature.py --sam /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-hi_uniquelyAligned.sam --info /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/87-hi.molbc.2.tsv --infoname molecular_bc --flagname M --out /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-hi_uniquelyAligned_molbc.sam

### 212-lo
python /g/steinmetz/project/mTEC_Seq/src/addFeature.py --sam /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/212-lo_uniquelyAligned.sam --info /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/212-lo.molbc.2.tsv --infoname molecular_bc --flagname M --out /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/212-lo_uniquelyAligned_molbc.sam


### 212-hi
python /g/steinmetz/project/mTEC_Seq/src/addFeature.py --sam /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/212-hi_uniquelyAligned.sam --info /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/212-hi.molbc.2.tsv --infoname molecular_bc --flagname M --out /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/212-hi_uniquelyAligned_molbc.sam

### 214-lo
python /g/steinmetz/project/mTEC_Seq/src/addFeature.py --sam /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/214-lo_uniquelyAligned.sam --info /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/214-lo.molbc.2.tsv --infoname molecular_bc --flagname M --out /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/214-lo_uniquelyAligned_molbc.sam

### 214-hi
python /g/steinmetz/project/mTEC_Seq/src/addFeature.py --sam /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/214-hi_uniquelyAligned.sam --info /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/sorting/214-hi.molbc.2.tsv --infoname molecular_bc --flagname M --out /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/214-hi_uniquelyAligned_molbc.sam



############## 
############## BUILD A BARCODE COLLAPSED SAM FILE (FULL READ)
############## GET TABLES ON UNIQUE MOLECULES AND UNIQUE BARCODES
############## 

# Copy the original version of script from
# /g/steinmetz/project/TSES/src/5Seq_rename-and-collapse.py to /g/steinmetz/project/mTEC_Seq/src/5Seq_rename-and-collapse.py
# Note: Minor edits on collapsing. 


### 87-lo
nohup python /g/steinmetz/project/mTEC_Seq/src/5Seq_rename-and-collapse.py --samfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-lo_uniquelyAligned_molbc.sam --out-collapsed /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/87-lo_uniquelyAligned_molbc_collapsed.sam --out-molecule-count /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/87-lo_collapsed_countsPerMol.tsv --out-position-count /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/87-lo_collapsed_countsPerPos.tsv --out-info /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/87-lo_collapsed.info > src/nohup.rename-and-collapse.87-lo.out &
[1] 14215
#--out-renamed /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-lo_uniquelyAligned_molbc_renamed.sam


### 87-hi
nohup python /g/steinmetz/project/mTEC_Seq/src/5Seq_rename-and-collapse.py --samfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/87-hi_uniquelyAligned_molbc.sam --out-collapsed /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/87-hi_uniquelyAligned_molbc_collapsed.sam --out-molecule-count /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/87-hi_collapsed_countsPerMol.tsv --out-position-count /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/87-hi_collapsed_countsPerPos.tsv --out-info /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/87-hi_collapsed.info > src/nohup.rename-and-collapse.87-hi.out &
[2] 14833


### 212-lo
nohup python /g/steinmetz/project/mTEC_Seq/src/5Seq_rename-and-collapse.py --samfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/212-lo_uniquelyAligned_molbc.sam --out-collapsed /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/212-lo_uniquelyAligned_molbc_collapsed.sam --out-molecule-count /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/212-lo_collapsed_countsPerMol.tsv --out-position-count /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/212-lo_collapsed_countsPerPos.tsv --out-info /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/212-lo_collapsed.info > src/nohup.rename-and-collapse.212-lo.out &
[3] 15233


### 212-hi
nohup python /g/steinmetz/project/mTEC_Seq/src/5Seq_rename-and-collapse.py --samfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/212-hi_uniquelyAligned_molbc.sam --out-collapsed /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/212-hi_uniquelyAligned_molbc_collapsed.sam --out-molecule-count /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/212-hi_collapsed_countsPerMol.tsv --out-position-count /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/212-hi_collapsed_countsPerPos.tsv --out-info /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/212-hi_collapsed.info > src/nohup.rename-and-collapse.212-hi.out &
[4] 15474


### 214-lo
nohup python /g/steinmetz/project/mTEC_Seq/src/5Seq_rename-and-collapse.py --samfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/214-lo_uniquelyAligned_molbc.sam --out-collapsed /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/214-lo_uniquelyAligned_molbc_collapsed.sam --out-molecule-count /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/214-lo_collapsed_countsPerMol.tsv --out-position-count /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/214-lo_collapsed_countsPerPos.tsv --out-info /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/214-lo_collapsed.info > src/nohup.rename-and-collapse.214-lo.out &
[5] 15636


### 214-hi
nohup python /g/steinmetz/project/mTEC_Seq/src/5Seq_rename-and-collapse.py --samfile /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/214-hi_uniquelyAligned_molbc.sam --out-collapsed /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/214-hi_uniquelyAligned_molbc_collapsed.sam --out-molecule-count /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/214-hi_collapsed_countsPerMol.tsv --out-position-count /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/214-hi_collapsed_countsPerPos.tsv --out-info /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/214-hi_collapsed.info > src/nohup.rename-and-collapse.214-hi.out &
[6] 15786


## The libraries have a lot of PCR amplification.


############## 
############## MAKE FILES FOR VISUALIZATION:
############## COLLAPSE READS TO THEIR MOST 5' POSITION; PREPARE SAM AND BEDGRAPGH
############## 

# Copy script /g/steinmetz/project/3Stability/src/divideByStrand.py 
# to /g/steinmetz/project/mTEC_Seq/src/

# Split by strand:
for i in collapsing/*_uniquelyAligned_molbc_collapsed.sam; do echo $i; python /g/steinmetz/project/mTEC_Seq/src/divideByStrand.py $i collapsing/`basename $i .sam`_plus.sam collapsing/`basename $i .sam`_minus.sam; done

# Sam to bam:
for i in collapsing/*_uniquelyAligned_molbc_collapsed_*.sam; do echo $i; samtools view -bS -o collapsing/`basename $i .sam`.bam $i; done

# Sort bam file
for i in collapsing/*_uniquelyAligned_molbc_collapsed_*us.bam; do echo $i; samtools sort $i collapsing/`basename $i .bam`.sorted; done

# Remove the sam
# rm collapsing/*_uniquelyAligned_molbc_collapsed*us.sam

# Index files (just for later browsing)
for i in collapsing/*_uniquelyAligned_molbc_collapsed_*us.sorted.bam; do echo $i; samtools index $i; done

### Collapse
for i in collapsing/*_uniquelyAligned_molbc_collapsed_*us.sorted.bam; do echo $i; genomeCoverageBed -5 -dz -ibam $i -g /g/steinmetz/genome/Homo_sapiens/37.68/fasta/DNA/Homo_sapiens.GRCh37.68.dna.chromosomes.withIVTs-chromLens.tsv > collapsing/`basename $i .sorted.bam`_5p.dz-collapsed.bedgraph; done



#####
##### Merge plus and minus strand
#####

# Copy script from g/steinmetz/project/PROMPT-seq/src/merge-bedgraphs.sh
# to /g/steinmetz/project/mTEC_Seq/src/


### 87-lo
sh /g/steinmetz/project/mTEC_Seq/src/merge-bedgraphs.sh /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing 87-lo_uniquelyAligned_molbc_collapsed_plus_5p.dz-collapsed.bedgraph 87-lo_uniquelyAligned_molbc_collapsed_minus_5p.dz-collapsed.bedgraph 87-lo_molbc-collapsed_5-collapsed.bedgraph

### 87-hi
sh /g/steinmetz/project/mTEC_Seq/src/merge-bedgraphs.sh /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing 87-hi_uniquelyAligned_molbc_collapsed_plus_5p.dz-collapsed.bedgraph 87-hi_uniquelyAligned_molbc_collapsed_minus_5p.dz-collapsed.bedgraph 87-hi_molbc-collapsed_5-collapsed.bedgraph

### 212-lo
sh /g/steinmetz/project/mTEC_Seq/src/merge-bedgraphs.sh /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing 212-lo_uniquelyAligned_molbc_collapsed_plus_5p.dz-collapsed.bedgraph 212-lo_uniquelyAligned_molbc_collapsed_minus_5p.dz-collapsed.bedgraph 212-lo_molbc-collapsed_5-collapsed.bedgraph

### 212-hi
sh /g/steinmetz/project/mTEC_Seq/src/merge-bedgraphs.sh /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing 212-hi_uniquelyAligned_molbc_collapsed_plus_5p.dz-collapsed.bedgraph 212-hi_uniquelyAligned_molbc_collapsed_minus_5p.dz-collapsed.bedgraph 212-hi_molbc-collapsed_5-collapsed.bedgraph

### 214-lo
sh /g/steinmetz/project/mTEC_Seq/src/merge-bedgraphs.sh /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing 214-lo_uniquelyAligned_molbc_collapsed_plus_5p.dz-collapsed.bedgraph 214-lo_uniquelyAligned_molbc_collapsed_minus_5p.dz-collapsed.bedgraph 214-lo_molbc-collapsed_5-collapsed.bedgraph

### 214-hi
sh /g/steinmetz/project/mTEC_Seq/src/merge-bedgraphs.sh /g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing 214-hi_uniquelyAligned_molbc_collapsed_plus_5p.dz-collapsed.bedgraph 214-hi_uniquelyAligned_molbc_collapsed_minus_5p.dz-collapsed.bedgraph 214-hi_molbc-collapsed_5-collapsed.bedgraph


## Format to proper bedgraph:

cat /g/steinmetz/project/mTEC_Seq/src/toBedgraph.R | R --slave --args dir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/' pattern='_5-collapsed.bedgraph'




############## 
############## GETTING YIELD NUMBERS
############## 

cat /g/steinmetz/project/mTEC_Seq/src/getPerformace.R | R --slave --args
indirAlignment='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/alignment/'
indirCollapsing ='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/'
inAlignmentPattern='uniquelyAligned.info'
inAmplificationPattern='_countsPerMol.tsv'
outDir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/performance/'
sortingOut='sorting-and-alignment.tsv'
amplificationOut='molecules-and-positions.tsv'
summaryOut='summary.tsv'



####
#### PLOTTING "EFFICIENCY"
####

cat /g/steinmetz/project/mTEC_Seq/src/plotEfficiency.R | R --slave --args indir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/performance/' infile='summary.tsv' outPdf='mTEC_5Seq_yield.pdf'


####
#### PLOTTING PCR AMPLIFICATION
####

cat /g/steinmetz/project/mTEC_Seq/src/plotPCRamplification.R | R --slave --args inDir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/' outDir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/performance/' inPattern='_countsPerMol.tsv' outPdf='PCR_amplification.pdf'

(Note: Don't mind the warnings)



############## 
############## SIMULATION: HOW CLOSE TO REACHING "SATURATION" (All 5' present)
############## 

cat /g/steinmetz/project/mTEC_Seq/src/estimateSaturation.R | R --slave --args inDir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/' outDir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/performance/' inPattern='_countsPerMol.tsv' outPdf='Simulation_coverage-to-identified-5-ends_all-data.pdf'






############## 
############## NUMBERS ON IVTS
############## 

cat /g/steinmetz/project/mTEC_Seq/src/getIVTaccuracy.R | R --slave --args infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv' indir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/' pattern='countsPerMol.tsv' outPlot='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/performance/5Seq-mapping-of-IVT-TSS.pdf'

## Strange: Mapping to IVT start sites looks like there is a lot of degradation (large fraction of counts not at TSS).





############## 
############## BUILD A COMMON TABLE OF 5' POSITIONS (SINGLE NT RESOLUTION)
############## 

cat /g/steinmetz/project/mTEC_Seq/src/make5positionTable.R | R --slave --args
indir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/'
pattern='_collapsed_countsPerMol.tsv'
outdir='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/'
outfile='collapsed-counts-table_genomic-chr_single-nt.tsv'





############## 
############## MAKE VENNS TO DISPLAY SAMPLE OVERLAP
############## 

# Note: plot5positionVenn.R works for up to 5 samples.

#### PER INDIVIDUAL

### 212 

## Cutoff 1 (all data)
cat /g/steinmetz/project/mTEC_Seq/src/plot5positionVenn.R | R --slave --args infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv' outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/plots/venn_5pos-overlap_indiv-212_cutoff1_single-nt_all-counts.tiff' samples='212_hi,212_lo' cutoff=1

## Cutoff 5
cat /g/steinmetz/project/mTEC_Seq/src/plot5positionVenn.R | R --slave --args
infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv'
outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/plots/venn_5pos-overlap_indiv-212_cutoff5_single-nt_all-counts.tiff'
samples='212_hi,212_lo'
cutoff=5

### 214

## Cutoff 1 (all data)
cat /g/steinmetz/project/mTEC_Seq/src/plot5positionVenn.R | R --slave --args infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv' outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/plots/venn_5pos-overlap_indiv-214_cutoff1_single-nt_all-counts.tiff' samples='214_hi,214_lo' cutoff=1

## Cutoff 5
cat /g/steinmetz/project/mTEC_Seq/src/plot5positionVenn.R | R --slave --args infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv' outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/plots/venn_5pos-overlap_indiv-214_cutoff5_single-nt_all-counts.tiff' samples='214_hi,214_lo' cutoff=5


### 87 

## Cutoff 1 (all data)
cat /g/steinmetz/project/mTEC_Seq/src/plot5positionVenn.R | R --slave --args infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv' outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/plots/venn_5pos-overlap_indiv-87_cutoff1_single-nt_all-counts.tiff' samples='87_hi,87_lo' cutoff=1

## Cutoff 5
cat /g/steinmetz/project/mTEC_Seq/src/plot5positionVenn.R | R --slave --args infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv' outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/plots/venn_5pos-overlap_indiv-87_cutoff5_single-nt_all-counts.tiff' samples='87_hi,87_lo' cutoff=5



#### AIRE LO

## Cutoff 1 (all data)
cat /g/steinmetz/project/mTEC_Seq/src/plot5positionVenn.R | R --slave --args infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv' outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/plots/venn_5pos-overlap_aire-lo_cutoff1_single-nt_all-counts.tiff' samples='212_lo,214_lo,87_lo' cutoff=1

## Cutoff 5
cat /g/steinmetz/project/mTEC_Seq/src/plot5positionVenn.R | R --slave --args infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv' outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/plots/venn_5pos-overlap_aire-lo_cutoff5_single-nt_all-counts.tiff' samples='212_lo,214_lo,87_lo' cutoff=5


#### AIRE HI

## Cutoff 1 (all data)
cat /g/steinmetz/project/mTEC_Seq/src/plot5positionVenn.R | R --slave --args infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv' outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/plots/venn_5pos-overlap_aire-hi_cutoff1_single-nt_all-counts.tiff' samples='212_hi,214_hi,87_hi' cutoff=1

## Cutoff 5
cat /g/steinmetz/project/mTEC_Seq/src/plot5positionVenn.R | R --slave --args infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv' outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/plots/venn_5pos-overlap_aire-hi_cutoff5_single-nt_all-counts.tiff' samples='212_hi,214_hi,87_hi' cutoff=5





############
############ COUNT CORRELATION PLOTS (SINGLE NT RESOLUTION)
############

## Without colors
cat /g/steinmetz/project/mTEC_Seq/src/plot5positionCountCorrelation.R | R --slave --args infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv' outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/plots/5pos-count-correlation.tiff' samples='212_hi,212_lo,214_hi,214_lo,87_hi,87_lo'


## LSD package style.

cat /g/steinmetz/project/mTEC_Seq/src/plot5positionCountCorrelation_LSD.R | R --slave --args infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv' outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/plots/5pos-count-correlation_LSD.tiff' samples='212_hi,212_lo,214_hi,214_lo,87_hi,87_lo'





############
############ TSS PLOTS
############

### GET TSS ANNOTATION FROM GTF FILE (IRANGES OBJECT)

cat /g/steinmetz/project/mTEC_Seq/src/getTSSfromGTF.R | R --slave --args
gtf='/g/steinmetz/genome/Homo_sapiens/37.68/annotation/gtf/Homo_sapiens.GRCh37.68.chrOnly.gtf'
outfile='/g/steinmetz/project/mTEC_Seq/annotation/TSS_IRanges.rda'


### PLOT 5SEQ SIGNAL AT TSSs

cat /g/steinmetz/project/mTEC_Seq/src/plot5positionAtTSS.R | R --slave --args infile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv' outMean='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/plots/Tx-TSS_TPM_w200.pdf' outHeatmap='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/plots/Tx-TSS_TPM_w200_heatmap.tiff' TSSIR='/g/steinmetz/project/mTEC_Seq/annotation/TSS_IRanges.rda' samples='212_hi,212_lo,214_hi,214_lo,87_hi,87_lo' window=200 dropTop=0.90




############
############ ASSIGNMENT OF 5' READS
############

####### Extract annotation features from GTF.

cat /g/steinmetz/project/mTEC_Seq/src/prepareAssignmentObjects.R | R --slave --args
gtf='/g/steinmetz/genome/Homo_sapiens/37.68/annotation/gtf/Homo_sapiens.GRCh37.68.chrOnly.gtf'
outfile='/g/steinmetz/project/mTEC_Seq/annotation/assignment_IRanges.rda'


####### Assign and plot pies.

cat /g/steinmetz/project/mTEC_Seq/src/makeAssignmentPies.R | R --slave --args data='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv' assignment='/g/steinmetz/project/mTEC_Seq/annotation/assignment_IRanges.rda' extension=100 samples='212_hi,212_lo,214_hi,214_lo,87_hi,87_lo' outPdfPosition='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/plots/5Seq-assignment-positions_TSS-100.pdf' outPdfCount='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/plots/5Seq-assignment-counts_TSS-100.pdf'





############## 
############## BROWSING THE DATA
############## 

#### Open IGV:

java -Xmx8000m -jar /g/steinmetz/jaerveli/software/IGV/IGV_2.3.7/igv.jar





# ----------------------------------------------------------> # From here on 28 March 2014.

############## 
############## PREPARING A GENE ANNOTATION OBJECT
############## 

cat /g/steinmetz/project/mTEC_Seq/src/prepareGeneAssignmentObject.R | R --slave --args gtf='/g/steinmetz/genome/Homo_sapiens/37.68/annotation/gtf/Homo_sapiens.GRCh37.68.chrOnly.gtf' outfile='/g/steinmetz/project/mTEC_Seq/annotation/genes-by-exons-plus-200-TSS-extension_IRanges.rda'



############## 
############## ASSIGNMENT TO GENES
############## 

cat /g/steinmetz/project/mTEC_Seq/src/assing2genes.R | R --slave --args data='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/collapsing/collapsed-counts-table_genomic-chr_single-nt.tsv' assignment='/g/steinmetz/project/mTEC_Seq/annotation/genes-by-exons-plus-200-TSS-extension_IRanges.rda' samples='212_hi,212_lo,214_hi,214_lo,87_hi,87_lo' outfile='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/expression/5Seq-at-genes.tsv'



############## 
############## DIFFERENTIAL GENE EXPRESSION
############## 

#### Install DESeq2.

cd /g/steinmetz/jaerveli/bin/R
wget http://www.bioconductor.org/packages/2.12/bioc/src/contrib/DESeq2_1.0.19.tar.gz
tar -zxvf DESeq2_1.0.19.tar.gz
# in R:
install.packages('/g/steinmetz/jaerveli/bin/R/DESeq2', lib = '/g/steinmetz/jaerveli/bin/R/DESeq2', repos= NULL, type="source")



#### RUN DESEQ ON GENES AND MAKE A FEW PLOTS

# Note: For DESeq2 to work, need a new R version (non-default R)
cat /g/steinmetz/project/mTEC_Seq/src/getDifferentiallyExpressedGenes.R | /g/software/bin/R-3.0.2 --slave --args geneTable='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/expression/5Seq-at-genes.tsv' outSig='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/expression/differentialExpression/genes-high-to-low-significant.tsv' outAll='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/expression/differentialExpression/genes-high-to-low-all.tsv' samples='212_hi,212_lo,214_hi,214_lo,87_hi,87_lo' groupIDs='hi,lo,hi,lo,hi,lo' outPlot='/g/steinmetz/project/mTEC_Seq/run/5TSeq_mTec1_2014-03-18_C3V0VACXX/assignment/expression/differentialExpression/genes-high-to-low-MA.png' sigCO=0.05



# ----------------------------------------------------------> # From here on DD Month YYYY.

## TO BE CONTINUED WITH ...


############## 
############## DIFFERENTIAL PROMOTER USAGE
############## 



############## 
############## ENRICHMENT ANALYSES OF DIFFERENTIALLY EXPRESSED GENES AND ISOFORMS
############## 



############## 
############## ANALYSIS OF TRANSCRIPTIONAL NOISE (AIRE- VS AIRE+)
############## 

## Promoter spread

## Promoter use (#)

## Intergenic expression



############## 
############## ANALYSIS OF AIRE CHIP-seq (& PolII ChIP-seq)??
############## 




############## 
############## PROCESSING CAGE DATA FROM FANTOM5 (TISSUE SPECIFIC CAGE)
############## 



############## 
############## FINDING MISSING EPITOPES (COMPARISON TO FANTOM5)
############## 





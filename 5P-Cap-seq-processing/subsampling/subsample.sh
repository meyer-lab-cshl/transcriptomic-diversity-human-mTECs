#!/bin/bash


function usage()
{
  cat <<EOF
  Usage: $0
          [ -c subsampling_count ]
          [ -s sample ]
          [ -r reads ]
          [ -t threads ]
          [ -i input_file ]
          [ -o output_file ]"

EOF
  exit 1;
}

while getopts c:s:r:t:i:o: opt; do
  case ${opt} in
  c) counts=${OPTARG};;
  s) sample=${OPTARG};;
  r) reads=${OPTARG};;
  t) threads=${OPTARG};;
  i) input=${OPTARG};;
  o) output=${OPTARG};;
  *) usage;;
  esac
done

#FACTOR=$(samtools idxstats $1 | cut -f3 | awk -v COUNT=$2 'BEGIN {total=0} {total += $1} END {print COUNT/total}')

echo "awk -v SAMPLE=$sample -v READS=$reads '(\$1 == SAMPLE) {print READS/\$2}' $counts"
FACTOR=$(awk -v SAMPLE=$sample -v READS=$reads '($1 == SAMPLE) {print READS/$2}' $counts)

#if [[ $FACTOR > 1 ]]
#  then
#  echo '[ERROR]: Requested number of reads exceeds total read count in' $1 '-- exiting' && exit 1
#fi

echo "sambamba view -s $FACTOR -t $threads --subsampling-seed=$reads -f bam -l 5 -o $output $input"
sambamba view -s $FACTOR -t $threads --subsampling-seed=$reads -f bam -l 5 -o $output $input


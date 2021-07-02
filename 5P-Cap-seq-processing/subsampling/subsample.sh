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

echo "awk -v SAMPLE=$sample -v READS=$reads '(\$1 == SAMPLE && \$2 == READS) {print $4}' $counts"
FACTOR=$(awk -v SAMPLE=$sample -v READS=$reads '($1 == SAMPLE && $2 == READS) {print $4}' $counts)

echo "sambamba view -s $FACTOR -t $threads --subsampling-seed=$reads -f bam -l 5 -o $output $input"
sambamba view -s $FACTOR -t $threads --subsampling-seed=$reads -f bam -l 5 -o $output $input


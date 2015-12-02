#!/bin/bash

# Sort a fastq file

awk 'ORS=(NR%4)?"}":"\n"' $1 > $1_tmp
rm -rf $1
sort -k 1 -t '}' $1_tmp | tr '}' '\n' > $1
rm -rf $1_tmp

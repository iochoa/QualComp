#!/bin/bash

# Parsing the arguments
numClusters=0;
numRates=0;
path=0;

while getopts p:c:r:h opt
do
  case "$opt" in
    c) clusters[numClusters]=$OPTARG
       (( numClusters = numClusters + 1));;
    r) rates[numRates]=$OPTARG
       (( numRates = numRates + 1 ));;
    p) path=$OPTARG;;
    h) echo " "
       echo "Options:"
       echo "-p <prefix of files to be decompressed>"
       echo "-r <rate in numBits/read>"
       echo "-c <number of clusters>"
       echo " "
       echo "Usage:"
       echo "Option -p: Needs to be used exactly one time"
       echo "Option -r: Needs to be used at least one time"
       echo "Option -c: (optional) It can be used multiple times. When not used, no cluster operation is assumed"
       echo "For multiple entries of -r or/and -c the program will generate several files for the differnt combinations of rate and number of clusters"
       echo " "
       echo "Example:"
       echo "./runDecompress.sh -p <prefix> -c 1 -c 4 -r 0 -r 100"
       echo "This will decompress the file with rate 0 bits/read and 100 bits/read (one cluster) and with 4 clusters and rate 0 bits/read and 100 bits/read" 
       echo " "       
       exit 1;;
  esac
done

# Checking if all the arguments are present
if [ "$numRates" == 0 ] || [ "$path" == 0 ] ; then
  echo "Error: You need to use options -p and -r at least once"
  echo "Use -h for help"
  exit 1
fi

if [ "$numClusters" == 0 ] ; then
  numClusters= 1;
  clusters[0]=1;
fi

# Print the arguments
echo "##########################"
echo "Prefix is" $path

for i in "${clusters[@]}"
do
  for rate in "${rates[@]}"
  do
    echo "##########################"
    echo "Computing MSE for" $i "clusters and rate" $rate
    ./computeMSE ${path}.fastq ${path}_${i}_${rate}.fastq `sed -n 1p ${path}_info` `sed -n 2p ${path}_info`
  done
done

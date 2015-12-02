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
echo "##########################"
echo "Number of clusters is" $numClusters
(( numClusters = numClusters - 1 ))
for i in $(seq 0 $numClusters)
do
  echo ${clusters[$i]}
done

echo "##########################"
echo "Number of rates is" $numRates
(( numRates = numRates - 1 ))
for i in $(seq 0 $numRates)
do
  echo ${rates[$i]}
done

# Start decompression 

# Iterate for each cluster and rate and decompress
for i in "${clusters[@]}"
do
  for j in $(seq 1 $i)
  do
    # Compute SVD
    ./computeSvd ${path}_${i}_${j}_stats `sed -n 1p ${path}_${i}_${j}_info` ${path}_${i}_${j}_svd 
  done

  for rate in "${rates[@]}"
  do
    echo "##########################"
    echo "Decompressing file with" $i "clusters and rate" $rate
    for j in $(seq 1 $i)
    do
      # Compute Rho
      ./computeOptimalRho ${path}_${i}_${j}_svd `sed -n 1p ${path}_${i}_${j}_info` $rate ${path}_${i}_${j}_${rate}_rho  

      # Decompress File
      ./decompressQval ${path}_${i}_${j}_${rate}_bin ${path}_${i}_${j}_stats ${path}_${i}_${j}_svd ${path}_${i}_${j}_${rate}_rho `sed -n 2p ${path}_${i}_${j}_info` `sed -n 1p ${path}_${i}_${j}_info` ${path}_${i}_${j}_${rate}_out `sed -n 3p ${path}_${i}_${j}_info` `sed -n 4p ${path}_${i}_${j}_info` 
      
      rm -rf ${path}_${i}_${j}_${rate}_rho

      # Create Fastq
      ./createFastq ${path}_noQvals_${i}_${j} ${path}_${i}_${j}_${rate}_out ${path}_${i}_${j}_${rate}.fastq
      rm -rf ${path}_${i}_${j}_${rate}_out
    done
    touch ${path}_${i}_${rate}.fastq_tmp
    for j in $(seq 1 $i)
    do
      cat ${path}_${i}_${j}_${rate}.fastq >> ${path}_${i}_${rate}.fastq_tmp
      rm -rf ${path}_${i}_${j}_${rate}.fastq
    done

    # Sort file
    ./sortFastq.sh ${path}_${i}_${rate}.fastq_tmp
    
    # Set Qvalue of N to minQV
    ./NtoMinQV ${path}_${i}_${rate}.fastq_tmp ${path}_${i}_${rate}.fastq `sed -n 1p ${path}_info` `sed -n 3p ${path}_info`
    rm -rf ${path}_${i}_${rate}.fastq_tmp

  done

  for j in $(seq 1 $i)
  do
    rm -rf ${path}_${i}_${j}_svd
  done

done

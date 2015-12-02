#!/bin/bash

# 1: Fastq file to be compressed
# 2: Length of the reads
# 3: Number of reads
# 4: File to store mean and covariance
# 5: File to store SVD information
# 6: Rate
# 7: File to store optimal Rho
# 8: Binary file where compressed file is stored
# 9: File with the uncompressed Quality Values

# Parsing the arguments
numClusters=0;
numRates=0;
inputFile=0;

while getopts nc:nr:i:c:r:h opt
do
  case "$opt" in
    c) clusters[numClusters]=$OPTARG
       (( numClusters = numClusters + 1));;
    r) rates[numRates]=$OPTARG
       (( numRates = numRates + 1 ));;
    i) inputFile=$OPTARG;;
    h) echo " "
       echo "Options:"
       echo "-i <input fastq file>"
       echo "-r <rate in numBits/read>"
       echo "-c <number of clusters>"
       echo " "
       echo "Usage:"
       echo "Option -i: Needs to be used exactly one time"
       echo "Option -r: Needs to be used at least one time"
       echo "Option -c: (optional) It can be used multiple times. When not used, no cluster operation is performed"
       echo "For multiple entries of -r or/and -c the program will compress the <fastq file> several times for the differnt combinations of rate and number of clusters"
       echo " "
       echo "Example:"
       echo "./runCompress.sh -i <fastq file> -c 1 -c 4 -r 0 -r 100"
       echo "This will compress the <fastq file> with rate 0 bits/read and 100 bits/read (one cluster) and it will also split the <fastq file> in 4 clusters, and compress each of them with 0 bits/read and 100 bits/read" 
       echo " "       
       exit 1;;
  esac
done

# Checking if all the arguments are present
if [ "$numRates" == 0 ] || [ "$inputFile" == 0 ] ; then
  echo "Error: You need to use options -i and -r at least once"
  echo "Use -h for help"
  exit 1
fi

if [ "$numClusters" == 0 ] ; then
  numClusters= 1;
  clusters[0]=1;
fi

# Print the arguments
echo "##########################"
echo "Fastq file is" $inputFile
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

# Extract filename, extension and path
filename=$(basename "$inputFile")
extension="${filename##*.}"
filename="${filename%.*}"
path="${inputFile%.*}"

# Start compression

# Separate file in QVs and no QVs
echo "##########################"
echo "Separating file in QVs and no QVs"
./separateQvals.sh $inputFile ${path}_Qvals ${path}_noQvals

# Compute numReads, szReads and minQV
echo "##########################"
echo "Computing numReads, szReads and minQV"
./computeNumSzReads ${path}_Qvals ${path}_info

# Iterate for each cluster and rate and compress
for i in "${clusters[@]}"
do
  echo "##########################"
  echo "Splitting data in" $i "clusters"
  # Run cluster
  if [ "$i" == 1 ] ; then
    cp ${path}_Qvals ${path}_Qvals_1
    cp ${path}_noQvals ${path}_noQvals_1
  else
    echo "$i is:" $i
    ./cluster ${path}_Qvals $i `sed -n 3p ${path}_info` `sed -n 4p ${path}_info` `sed -n 1p ${path}_info` `sed -n 2p ${path}_info` ${path}_noQvals
  fi

  for j in $(seq 1 $i)
  do
    mv ${path}_noQvals_${j} ${path}_noQvals_${i}_${j}
    mv ${path}_Qvals_${j} ${path}_Qvals_${i}_${j}
    
    # Compute mean and covariance
    ./computeNumSzReads ${path}_Qvals_${i}_${j} ${path}_${i}_${j}_info
    ./generateStats ${path}_Qvals_${i}_${j} `sed -n 1p ${path}_${i}_${j}_info` `sed -n 2p ${path}_${i}_${j}_info` ${path}_${i}_${j}_stats

    # Compute SVD
    ./computeSvd ${path}_${i}_${j}_stats `sed -n 1p ${path}_${i}_${j}_info` ${path}_${i}_${j}_svd 
  done

  for rate in "${rates[@]}"
  do
    echo "##########################"
    echo "Compressing file with" $i "clusters and rate" $rate
    for j in $(seq 1 $i)
    do
      # Compute Rho
      ./computeOptimalRho ${path}_${i}_${j}_svd `sed -n 1p ${path}_${i}_${j}_info` $rate ${path}_${i}_${j}_${rate}_rho  

      # Compress File
      if [ "$rate" == 0 ] ; then
        touch ${path}_${i}_${j}_${rate}_bin
      else
        ./compressQval ${path}_Qvals_${i}_${j} ${path}_${i}_${j}_stats ${path}_${i}_${j}_svd ${path}_${i}_${j}_${rate}_rho `sed -n 2p ${path}_${i}_${j}_info` `sed -n 1p ${path}_${i}_${j}_info` ${path}_${i}_${j}_${rate}_bin 
      fi

       rm -rf ${path}_${i}_${j}_${rate}_rho
    done
  done
  
  for j in $(seq 1 $i)
  do
    rm -rf ${path}_${i}_${j}_svd
    rm -rf ${path}_Qvals_${i}_${j}
  done

done

rm -rf ${path}_Qvals
rm -rf ${path}_noQvals
# rm -rf ${path}_info


# echo "Number of bytes of files is"
#ls -lhtr $8

# Uncompress the file
#./decompressQval $8 $4 $5 $7 $3 $2 $9

# Compute MSE
#./computeMSE $1_qvals $9 $2 $3

# Example of how to run the program

# Download the data 'SRR032209.fastq' and place it in a folder named 'SRR032209' in the same directory as the lossy compression algorithm

# Compress it using 1 and 2 clusters, and with 0 and 36 bits/quality score sequence
./runCompress.sh -i SRR032209/SRR032209.fastq -c 1 -c 2 -r 0 -r 36

# This will create, in the directory 'SRR032209', the following files:

# Info of option 1 cluster (no cluster)
# SRR032209_1_1.info
# SRR032209_1_1.stats

# Compressed files for option 1 cluster and rate 0
# SRR032209_1_1_0.bin

# Compressed files for option 1 cluster and rate 36
# SRR032209_1_1_36.bin


# Info of option 2 clusters
# SRR032209_2_1.info
# SRR032209_2_2.info
# SRR032209_2_1.stats
# SRR032209_2_2.stats

# Compressed files for option 2 clusters and rate 0
# SRR032209_2_1_0.bin
# SRR032209_2_2_0.bin

# Compressed files for option 2 clusters and rate 36
# SRR032209_2_1_36.bin
# SRR032209_2_2_36.bin

# Extra information contained in the fastq file (all except quality scores)
# SRR032209_noQvals_1_1
# SRR032209_noQvals_2_1
# SRR032209_noQvals_2_2

# Create new fastq files with the reconstructed quality scores
./runDecompress.sh -p SRR032209/SRR032209 -c 1 -c 2 -r 0 -r 36

# This will generate the following files

# SRR032209_1_0.fastq
# SRR032209_1_36.fastq
# SRR032209_2_0.fastq
# SRR032209_2_36.fastq

# To compute the MSE
./sortFastq.sh SRR032209/SRR032209.fastq
./computeMSE.sh -p SRR032209/SRR032209 -c 1 -c 2 -r 0 -r 36





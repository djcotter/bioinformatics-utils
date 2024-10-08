#!/bin/bash

# This script takes in a list of SRA or ENA accessions and uses ffq 
# to write a list of paths to the fastq files to stdout.

# The script requires the ffq tool to be installed and in the PATH.

# get the list of accessions
accessions=$1
output=$2

# clear the output file
> $output

# loop over the accessions and get the fastq files
# and write the paths to the output file
# do this in parallel
cat $accessions | parallel --delay .5 -j 32 ffq --ftp {} | grep -Po '(?<=url": ")ftp://.*(?=")' >> $output


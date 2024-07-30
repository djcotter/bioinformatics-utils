#!/bin/bash

# This script takes in a list of SRA or ENA accessions and uses ffq 
# to write a list of paths to the fastq files to stdout.

# The script requires the ffq tool to be installed and in the PATH.

# get the list of accessions
accessions=$1

# loop over the accessions and get the fastq files
while read acc; do
    ffq --ftp $acc | grep -Po  '(?<=url": ")ftp://.*(?=")'
done < $accessions


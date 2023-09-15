#!/usr/bin/env bash
# -*- coding:  utf-8 -*-
# Author: Gaurav Sablok
# date: 2023-09-15
# MIT License
# a complete pangenome pipeline from the graph generation to the mapping of the 
# reads and estimation of the variants and visualization on the panchae. This uses
# a bacterial metagenomics approach and if you want a complete workflow for the plant
# look at the other respository
echo    "I am a bacterial metagenome assembler from the raw metagenomics reads to the 
            creation of the metagenomes, annotations and metagenomes pangraphs
             and visualization of the linear metagenomes in panache"
echo        "PLEASE REMOVE THE ECHO WORD BEFORE SUBMITTING ON THE CLUSTER"

read -p "type of machine that you are running this workflow on:" machine
read -p "do you want to follow the workflow": answer
read -p "please provide the path to the directory with reads:" reads
read -p "please provide the working directory path on the cluster:" workdir
read -p "do you want to clean the reads:" clean
read -p "provide the index name:" index
read -p "provide the mapper name:" mapper
if [[ $machine == "desktop" ]]
then
    username=$(users)
        echo "creating environment"
        echo "conda create -n metagenomics -y && conda install -n metagenomics python=3.11 -y"
        echo "conda activate metagenomics"
        echo "conda install blasr"
        echo "conda install blast2"
        echo "conda install canu"
        echo "conda install flye"
        echo "conda deactivate"
        echo "sudo apt-get install bowtie2"
        echo "sudo apt-get install spades"
	echo "finish creating the environment for the metagenomics analysis"
    break
fi
if [[ $machine == "server" ]] 
then 
     read -p "please provide the username:" username
     read -p "please provide the ssh_address:" ssh_address
     echo "ssh ${username}@ssh_address"
     break
fi
if [[ $answer == "yes" ]] && 
            [[ -d $reads ]] &&
               [[ -d $workdir ]]
then 
  mkdir ${workdir}/reads ${workdir}/metagenome_assembly ${workdir}/cleaned_reads ${workdir}/sam_files ${workdir}/bam_files
  mv ${reads}/*.fastq ${workdir}/reads
fi
readpath=${workdir}/reads
dirname=${workdir}/cleaned_reads
if [[ -d "${readpath}" ]] &&
        [[ $clean == "yes" ]] &&
            [[ $mapper == "yes" ]] &&
                [[ $index == "$index" ]]
then
    for i in ${readpath}/*.R1.fastq; do
        for j in ${readpath}/*.R2.fastq; do
        echo "Processing the cleaning of the reads file from the sequencing runs"
        done
        echo "fastp --in1 $i --out1 $i.clean.R1.fastq --in2 $j --out2 $j.clean.R2.fastq --threads $thread"
    done
fi
    echo "moving the cleaned reads to the respective folder"
        mv ${readpath}/*.cleaned.R1.fastq ${dirname}/cleaned_reads
        mv ${readpath}/*.cleaned.R2.fastq ${dirname}/cleaned_reads
    for i in "${dirname}"/*.cleaned.R1.fastq; do
        for j in "${dirname}"/*cleaned.R2.fastq; do
        echo "Processing the mapping of the reads file from the sequencing runs"
        done
        echo "bowtie2-build $index $index"
        echo "bowtie2 -t -x $index -p $thread --very-sensitive-local -1 $i -2 $j -S $index.sam --no-unal --al-conc $index.aligned.fastq"
    done

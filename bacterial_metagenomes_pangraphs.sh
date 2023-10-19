#!/usr/bin/env bash
# -*- coding:  utf-8 -*-
# Author: Gaurav Sablok
# date: 2023-09-15
# updated: 2023-10-19 with support for the remapping, snpcalling
# next update will be the support for the multi mappers, graph creation, pangenome construction
# evolutionary analysis, genome comparison
# python code will be intergated, plus ruby bindings and also a 
# ruby on rails application for submission to the server
# MIT License
# a complete pangenome pipeline from the graph generation to the mapping of the 
# reads and estimation of the variants and visualization on the panchae. This uses
# a bacterial metagenomics approach and if you want a complete workflow for the plant
# look at the other respository
echo    "I am a bacterial metagenome assembler from the raw metagenomics reads to the 
            creation of the metagenomes, annotations and metagenomes pangraphs
             and visualization of the linear metagenomes in panache"
echo        "PLEASE REMOVE THE ECHO WORD BEFORE SUBMITTING ON THE CLUSTER"

read -r -p "type of machine that you are running this workflow on:" machine
read -r -p "please provide the threads": thread
read -r -p "please provide the path to the directory with reads:" reads
read -r -p "please provide the working directory path on the cluster:" workdir
read -r -p "do you want to clean the reads:" clean
read -r -p "provide the mapper name:" mapper
read -r -p "provide the index name:" index
read -r -p "provide the answer for the assembly:" option
read -r -p "do you want to call the snps:" snps
read -r -p "please provide the read depth for the snps:" readdepth
read -r -p "please provide the depth for the snps:" filter

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
fi

if [[ $machine == "server" ]] 
then 
     read -r -p "please provide the username:" username
     read -r -p "please provide the ssh_address:" ssh_address
     read -r -p "please provide the threads": thread
     read -r -p "please provide the path to the scratch directory:" directory
     read -r -p "please provide the path to the directory with reads:" reads
     read -r -p "please provide the working directory path on the cluster:" workdir
     read -r -p "do you want to clean the reads:" clean
     read -r -p "provide the mapper name:" mapper
     read -r -p "provide the index name:" index
     read -r -p "provide the answer for the assembly:" option
     read -r -p "do you want to call the snps:" snps
     read -r -p "please provide the read depth for the snps:" readdepth
     read -r -p "please provide the depth for the snps:" filter
     echo "ssh ${username}@${ssh_address}"
     echo ""cd "${directory}"""
     mkdir  "${workdir}"/reads \
          "${workdir}"/metagenome_assembly \
          "${workdir}"/cleaned_reads 
          "${workdir}"/sam_files \
          "${workdir}"/aligned_reads \
          "${workdir}"/remapping \
          "${workdir}"/snpcalling \
          readpath="${workdir}"/reads
          dirname="${workdir}"/cleaned_reads
          aligned="${workdir}"/aligned_reads
          genomeassembly="${workdir}"/metagenome_assembly
          remapping="${workdir}"/remapping
          snpcalling="${workdir}"/snpcalling
          for i in "${readpath}"/*.R1.fastq
            do
                for j in "${readpath}"/*.R2.fastq
            do
                echo "Processing the cleaning of the reads file from the sequencing runs"
          done
            echo "fastp --in1 $i --out1 $i.clean.R1.fastq --in2 \
                                    $j --out2 $j.clean.R2.fastq --threads $thread"
           done
            echo "moving the cleaned reads to the respective folder"
            mv "${readpath}"/*.cleaned.R1.fastq "${dirname}"/cleaned_reads
            mv "${readpath}"/*.cleaned.R2.fastq "${dirname}"/cleaned_reads
            for i in "${dirname}"/*.cleaned.R1.fastq
                do
                    for j in "${dirname}"/*cleaned.R2.fastq
                do
                    echo "Processing the mapping of the reads file from the sequencing runs"
                done
                    echo "bowtie2-build $index $index"
                    echo "bowtie2 -t -x $index -p $thread --very-sensitive-local -1 $i -2 $j \
                                            -S $index.sam --no-unal --al-conc $index.aligned.fastq"
                    echo ""mv *.aligned.fastq "${aligned}"""
                done
            for i in "${aligned}"/*.aligned.R1.fastq
                do
                    for j in "${aligned}"/*.aligned.R2.fastq
                do
                    echo "Processing genome assembly"
                done
                    echo "spades.py -1 $i -2 $j --careful --threads $thread --tmp-dir \
                                         ${index}_tmp -k 45,55,65,75 -o ${dirname}/${index}_assembly"
                    echo ""mv *.fasta "${genomeassembly}"""
                done

      echo "moving the reads for the remapping"
      echo ""cp -r "${genomeassembly}/*fasta" "${remapping}"""
      echo ""cp -r "${dirname}/*.fastq" "${remapping}"""
      echo "bowtie2-build $index $index"
      echo "bowtie2 -t -x $index -p $thread --very-sensitive-local -1 $i -2 $j \
                                            -S $index.sam --no-unal --al-conc $index.aligned.fastq"
      for f in *.sam; do echo samtools -in "{$f}" -out "{$f%}".bam; done
      for f in *.bam; do echo "{$f}"; bamtools stats -in \
                                        "{$f}" --insert >> "{$f%}".insert; done 
      for f in *.bam; do echo "{$f}"; bamtools coverage \
                                        -in "{$f}" -out "{$f%}".coverage.txt; done
      for f in *.bam; do echo "{$f}"; bamtools count -in "{$f}" \
                                        > "{$f%}".count.read.alignment.txt; done
      touch insert_coverage_alignment.txt | paste "{$f}".insert_size.txt \
                                        "{$f}".coverage.txt "{$f}".count.read.alignment.txt
        echo "mv ${remapping}"/*.sam "${snpcalling}"
        echo ""mv "${genomeassembly}"/*.fasta "${snpcalling}"""                                
        echo "read cleaning, mapping of the raw reads, genome assembly, 
                        remapping finished and processing to the snps calculation"
        echo ""cd "${snpcalling}""
        fastagenome=$(*.fasta)
        samfile=$(*.sam)
        samtools view -bs "${samfile} -o "${samfile‰}".bam
        samtools sort "${samfile‰}".bam "${samfile‰}".sorted.bam
        bamtools index "${samfile‰}".sorted.bam
        samtools index faidx "${fastagenome}"
        samtools mpileup -g -f "${fastagenome}" > "${fastagenome%}".bcf
        bcftools view -bvcg "${fastagenome%}".bcf > "${fastagenome%}".variant.bcf
        bcftools view "${fastagenome%}".variant.bcf | vcfutils.pl varFilter \
                        -d "${filter}" -d "${readdepth}" > "${fastagenome%}".selected.vcf
        gatk AnnotateVcfWithBamDepth --input "${samfile‰}".sorted.bam \
                                    --reference "${fastagenome}" \
                                    --output "${fastagenome}".annotated.vcf \
                                                    --variant "${fastagenome%}".selected.vcf
fi

if [[ $machine == "desktop" ]] &&
          [[ -d "${readpath}" ]] &&
              [[ $clean == "yes" ]] &&
                [[ $mapper == "yes" ]] &&
                     [[ $index == "$index" ]]
then 
  mkdir "${workdir}"/reads \
        "${workdir}"/cleaned_reads 
        "${workdir}"/sam_files \
        "${workdir}"/aligned_reads \
        mv ${reads}/*.fastq ${workdir}/reads
        readpath="${workdir}"/reads
        dirname="${workdir}"/cleaned_reads
        aligned="${workdir}"/aligned_reads
        for i in "${readpath}"/*.R1.fastq 
        do
            for j in "${readpath}"/*.R2.fastq
            do
                echo "Processing the cleaning of the reads file from the sequencing runs"
            done
                echo "fastp --in1 $i --out1 $i.clean.R1.fastq --in2 \
                                            $j --out2 $j.clean.R2.fastq --threads $thread"
        done
            echo "moving the cleaned reads to the respective folder"
                mv "${readpath}"/*.cleaned.R1.fastq "${dirname}"/cleaned_reads
                mv "${readpath}"/*.cleaned.R2.fastq "${dirname}"/cleaned_reads
        for i in "${dirname}"/*.cleaned.R1.fastq
            do
                for j in "${dirname}"/*cleaned.R2.fastq
            do
                echo "Processing the mapping of the reads file from the sequencing runs"
            done
                echo "bowtie2-build $index $index"
                echo "bowtie2 -t -x $index -p $thread --very-sensitive-local -1 $i -2 $j \
                                            -S $index.sam --no-unal --al-conc $index.aligned.fastq"
                echo ""mv *.aligned.fastq "${aligned}"""
        done
elif [[ $machine == "desktop" ]] &&
            [[ -d "${readpath}" ]] &&
                [[ $clean == "yes" ]] &&
                    [[ $mapper == "yes" ]] &&
                        [[ $index == "$index" ]] &&
                                [[ $option == "yes" ]] && 
                                        [[ $snps == "yes" ]] 
then
    mkdir  "${workdir}"/reads \
           "${workdir}"/metagenome_assembly \
           "${workdir}"/cleaned_reads 
           "${workdir}"/sam_files \
           "${workdir}"/aligned_reads \
           readpath="${workdir}"/reads
           dirname="${workdir}"/cleaned_reads
           aligned="${workdir}"/aligned_reads
           genomeassembly="${workdir}"/metagenome_assembly
      for i in "${readpath}"/*.R1.fastq; do
        for j in "${readpath}"/*.R2.fastq; do
            echo "Processing the cleaning of the reads file from the sequencing runs"
        done
            echo "fastp --in1 $i --out1 $i.clean.R1.fastq --in2 \
                                    $j --out2 $j.clean.R2.fastq --threads $thread"
        done
        echo "moving the cleaned reads to the respective folder"
            mv "${readpath}"/*.cleaned.R1.fastq "${dirname}"/cleaned_reads
            mv "${readpath}"/*.cleaned.R2.fastq "${dirname}"/cleaned_reads
        for i in "${dirname}"/*.cleaned.R1.fastq
            do
                for j in "${dirname}"/*cleaned.R2.fastq
            do
            echo "Processing the mapping of the reads file from the sequencing runs"
            done
                echo "bowtie2-build $index $index"
                echo "bowtie2 -t -x $index -p $thread --very-sensitive-local -1 $i -2 $j \
                                            -S $index.sam --no-unal --al-conc $index.aligned.fastq"
                echo ""mv *.aligned.fastq "${aligned}"""
            done
        for i in "${aligned}"/*.aligned.R1.fastq
            do
                for j in "${aligned}"/*.aligned.R2.fastq
            do
                echo "Processing genome assembly"
            done
            echo ""spades.py -1 $i -2 $j --careful --threads $thread --tmp-dir \
                                         "${index}_tmp" -k 45,55,65,75 -o "${dirname}/${index}_assembly""
            echo ""mv *.fasta "${genomeassembly}""
            done

elif [[ $machine == "desktop" ]] &&
            [[ -d "${readpath}" ]] &&
                [[ $clean == "yes" ]] &&
                    [[ $mapper == "yes" ]] &&
                        [[ $index == "$index" ]] &&
                                [[ $option == "yes" ]] &&
                                   [[ $snps == "yes" ]]                                
then 
   mkdir  "${workdir}"/reads \
          "${workdir}"/metagenome_assembly \
          "${workdir}"/cleaned_reads 
          "${workdir}"/sam_files \
          "${workdir}"/aligned_reads \
          "${workdir}"/remapping \
          "${workdir}"/snpcalling \
          readpath="${workdir}"/reads
          dirname="${workdir}"/cleaned_reads
          aligned="${workdir}"/aligned_reads
          genomeassembly="${workdir}"/metagenome_assembly
          remapping="${workdir}"/remapping
          snpcalling="${workdir}"/snpcalling
          for i in "${readpath}"/*.R1.fastq
            do
                for j in "${readpath}"/*.R2.fastq
            do
                echo "Processing the cleaning of the reads file from the sequencing runs"
          done
            echo "fastp --in1 $i --out1 $i.clean.R1.fastq --in2 \
                                    $j --out2 $j.clean.R2.fastq --threads $thread"
           done
            echo "moving the cleaned reads to the respective folder"
            mv "${readpath}"/*.cleaned.R1.fastq "${dirname}"/cleaned_reads
            mv "${readpath}"/*.cleaned.R2.fastq "${dirname}"/cleaned_reads
            for i in "${dirname}"/*.cleaned.R1.fastq
                do
                    for j in "${dirname}"/*cleaned.R2.fastq
                do
                    echo "Processing the mapping of the reads file from the sequencing runs"
                done
                    echo "bowtie2-build $index $index"
                    echo "bowtie2 -t -x $index -p $thread --very-sensitive-local -1 $i -2 $j \
                                            -S $index.sam --no-unal --al-conc $index.aligned.fastq"
                    echo ""mv *.aligned.fastq "${aligned}"""
                done
            for i in "${aligned}"/*.aligned.R1.fastq
                do
                    for j in "${aligned}"/*.aligned.R2.fastq
                do
                    echo "Processing genome assembly"
                done
                    echo "spades.py -1 $i -2 $j --careful --threads $thread --tmp-dir \
                                         ${index}_tmp -k 45,55,65,75 -o ${dirname}/${index}_assembly"
                    echo ""mv *.fasta "${genomeassembly}"""
                done

      echo "moving the reads for the remapping"
      echo ""cp -r "${genomeassembly}/*fasta" "${remapping}"""
      echo ""cp -r "${dirname}/*.fastq" "${remapping}"""
      echo "bowtie2-build $index $index"
      echo "bowtie2 -t -x $index -p $thread --very-sensitive-local -1 $i -2 $j \
                                            -S $index.sam --no-unal --al-conc $index.aligned.fastq"
      for f in *.sam; do echo samtools -in "{$f}" -out "{$f%}".bam; done
      for f in *.bam; do echo "{$f}"; bamtools stats -in \
                                        "{$f}" --insert >> "{$f%}".insert; done 
      for f in *.bam; do echo "{$f}"; bamtools coverage \
                                        -in "{$f}" -out "{$f%}".coverage.txt; done
      for f in *.bam; do echo "{$f}"; bamtools count -in "{$f}" \
                                        > "{$f%}".count.read.alignment.txt; done
      touch insert_coverage_alignment.txt | paste "{$f}".insert_size.txt \
                                        "{$f}".coverage.txt "{$f}".count.read.alignment.txt
        echo "mv ${remapping}"/*.sam "${snpcalling}"
        echo ""mv "${genomeassembly}"/*.fasta "${snpcalling}"""                                
        echo "read cleaning, mapping of the raw reads, genome assembly, 
                        remapping finished and processing to the snps calculation"
        echo "cd ${snpcalling}"
        fastagenome=()
        for file in "$(pwd)"/*.fasta
          do 
            fastagenome+=("file")
          done
            if [[ -z "${fastagenome}" ]]
                then
                    echo "file not present"
            fi
        samfile=()
        for file in "$(pwd)"/*.sam
           do 
                samfile+=("file")
           done
            if [[ -z "${samfile}" ]]
                then
                    echo "file not present"
                fi
        echo "samtools view -bs "${samfile}" -o "${samfile‰}".bam"
        echo "samtools sort "${samfile‰}".bam "${samfile‰}".sorted.bam"
        echo "bamtools index "${samfile‰}".sorted.bam"
        echo "samtools index faidx "${fastagenome}""
        echo "samtools mpileup -g -f "${fastagenome}" > "${fastagenome%}".bcf"
        echo "bcftools view -bvcg "${fastagenome%}".bcf > "${fastagenome%}".variant.bcf"
        echo "bcftools view "${fastagenome%}".variant.bcf | vcfutils.pl varFilter \
                        -d "${filter}" -d "${readdepth}" > "${fastagenome%}".selected.vcf"
        echo "gatk AnnotateVcfWithBamDepth --input "${samfile‰}".sorted.bam \
                                    --reference "${fastagenome}" \
                                    --output "${fastagenome}".annotated.vcf \
                                                    --variant "${fastagenome%}".selected.vcf"
fi                                                                                                       

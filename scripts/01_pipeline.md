---
title: Read trimming, metagenome assembly, and read recruitment
---

This pipeline starts with raw metagenome reads, (1) trims the raw reads, (2) assembles metagenomes from the trimmed reads, and (3) maps trimmed reads to the metagenomes.

## Trim reads

Create a NovaSeq.fa file containing the adapter sequences:
```
>PrefixPE/1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>PrefixPE/2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```
Trim adapter sequences, low quality bases, and ultrashort reads with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic "Trimmomatic") version 0.38:

```shell
java –jar trimmomatic-0.38.jar PE R1.fastq.gz R2.fastq.gz
R1_p_trimmed.fastq.gz R1_u_trimmed.fastq.gz R2_p_trimmed.fastq.gz R2_u_trimmed.fastq.gz
ILLUMINACLIP:NovaSeq.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
- Input files:
  - R1.fastq.gz and R2.fastq.gz are the raw, zipped .fastq files.
  - NovaSeq.fa is the .fasta file containing the adapter sequences.
- Output files:
  - R1_p_trimmed.fastq.gz and R2_p_trimmed.fastq.gz are the trimmed, paired reads.
  - R1_u_trimmed.fastq.gz and R2_u_trimmed.fastq.gz are the trimmed, unpaired reads.

#### Count trimmed reads

Count the number of trimmed, paired reads:
```shell
echo $(zcat R1_p_trimmed.fastq.gz|wc -l)/4|bc
```

## Assemble metagenomes

Assemble metagenomes with [MEGAHIT](https://github.com/voutcn/megahit "MEGAHIT") version 1.0.6 from the trimmed, paired reads:

```shell
megahit -1 R1_p_trimmed.fastq.gz -2 R2_p_trimmed.fastq.gz
--k-list 23,43,63,83,103,123 -o OUTPUT_DIRECTORY/ --verbose
```
- Output files:
  - contigs.fa is the metagenome assembly contigs file.

## Recruit reads to metagenomes

First, build a reference for the assembly:

```shell
bbmap.sh ref=contigs.fa path=metagenome_ref

```

Align reads at a 70 % sequence identity threshold. In the below example, bottom sediment metagenome reads (of Lac Paula) are mapped to the top sediment metagenome assembly (from the same lake):

```shell
cd lacpaula-top_ref/

bbmap.sh maxlen=500 ambiguous=best mappedonly=t idtag=t pigz=t minid=0.70
covstats=lacpaula-bottom_recruit_to_top_assembly_90id_covstats.txt
idhist=lacpaula-bottom_recruit_to_top_assembly_90id_idhist.txt
in=lacpaula-bottom_R1_p_trimmed.fastq.gz
in2=lacpaula-bottom_R2_p_trimmed.fastq.gz
outm=lacpaula-bottom_recruit_to_top_assembly_90id_mapped.sam
```

Align reads at a 90 % sequence identity threshold.  Again in the below example, bottom sediment metagenome reads (of Lac Paula) are mapped to the top sediment metagenome assembly (from the same lake):

```shell
cd lacpaula-top_ref/

bbmap.sh maxlen=500 ambiguous=best mappedonly=t idtag=t pigz=t minid=0.90
covstats=lacpaula-bottom_recruit_to_top_assembly_90id_covstats.txt
idhist=lacpaula-bottom_recruit_to_top_assembly_90id_idhist.txt
in=lacpaula-bottom_R1_p_trimmed.fastq.gz
in2=lacpaula-bottom_R2_p_trimmed.fastq.gz
outm=lacpaula-bottom_recruit_to_top_assembly_90id_mapped.sam
```

#### Count number mapped reads

Convert .sam to .bam:
```
samtools view -S -b 90id_mapped.sam > 90id_mapped.bam
```

Sort .bam:
```
samtools sort 90id_mapped.bam -o 90id_mapped_sorted.bam
```

Index sorted .bam:
```
samtools index -b 90id_mapped_sorted.bam
```

Write idxstats.txt file with [samtools idxstats](http://www.htslib.org/doc/samtools-idxstats.html]):
```
samtools idxstats 90id_mapped_sorted.bam > 90id_mapped_idxstats.txt
```
- idxstats.txt file contains the columns:
  - Reference sequence
  - Contig length
  - Number of mapped read-segments
  - Number of unmapped reads-segments

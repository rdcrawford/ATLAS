# ATLAS:  
Alignment of Targeted Locus Amplification Sequencing 

  Author:  Ryan D. Crawford
  
  Contact: rcrawfo@umich.edu


## Introduction: ##

  ATLAS (Alignment of Targeted Locus Amplification Sequencing) is designed for
  analysis of TLA data used in high throughput mutagenesis experiments. This
  program derives functionality for analyzing BAM files from the classes
  implimented in the bamtools library. First the raw paired reads are aligned
  with BWA mem with the penalty for unsparring reads set to zero. The the bam
  file is then parsed for reads that map to the locus containing a barcode
  sequence. Reads mapping to this locus are then searched against a library of
  barcode sequences and if a match is found, the opposite mate is accessed if it
  is duplicately mapped. Duplicately mapped reads are then digested with the
  restriction enzyme used to create the library, rewritten to fastq format, and
  remapped to the reference genome. The resulting BAM file is then tagged with
  the index of the barcode sequence contained in the opposite read. The BAM
  file is then sorted and written.


## Dependencies: ##

  1) Bamtools C++ API
  2) BWA
  3) Samtools



## Input parameters: ##

  1) Reference genome (fasta format)
     - If the sequence of interest is knocked in, the native locus should be
       mask from the reference and the sequence of the knocked in construct
       should be appended to the reference genome.
  2) Mate Pair 1 (fastq format)
  3) Mate pair 2 (fastq format)
  4) Barcode sequences (fasta format)
  5) Restriction site sequence
  6) Name of the chromosome containing the anchor sequence (from the fasta
     header)
  7) Barcode start position (0 indexed)
  8) Barcode end position (0 indexed)


## Output: ##

  1) Bamfile containing only informative reads. The index of the barcode
      contained in the opposite mate pair is indicated under the tag "XB".
  2) Tab delimited file containing the barcode sequences, the number assigned to
     each barcode and written in the bamfile under the "XB" tag.
  3) Fastq file containing the processed read. Some of which will be
     restriction fragments from the



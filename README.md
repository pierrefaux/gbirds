# GBirdS

GBirdS is a pipeline for avian population genetics from GBS reads.

See Faux et al. Fast genomic analysis of aquatic bird populations from short single-end reads considering sex-related pitfalls.


The current version of the pipeline (v1.2) includes two files:
  (1) gbirds12.f90 : the main body of the pipeline, including demultiplexing, sort of reads into unique sequences per sample and for the whole population, quality checks, identification of sex-linked sequences, genotype calling and filtering;
  (2) FinalSteps.m: a Matlab/Octave function including downstream steps: G matrix computation, Fst computation and principal component analysis.
The script "gbirds12.f90" includes all its dependencies (subroutines and functions) and was compiled on a Linux 64-bit machine in the binary file "gbirds".

After copying the binary files and making it accessible to your search path, please follow these steps to run the pipeline:
(1) Make a new work directory

(2) Copy the (unzipped) FASTQ files (sequencing lanes) to that directory

(3) Copy the population map (see below) to that directory

(4) Write a file containing the list of FASTQ, optionnaly followed in a 2nd (space-delimited) column by the number of reads in     each corresponding FASTQ 

(5) Optional: write an option file (see below for options)

(6) Run the pipeline with this command:

    gbirds --map [population map file] --files [FASTQ list file]
    
Or, if you have an option file to pass:
    
    gbirds --map [population map file] --files [FASTQ list file] --options [options file]
    
    

POPULATION FILE FORMAT:
$1: Sample barcode (up to 10 bp)
$2: ID (up to 20 characters)
$3: Population ID (integer, from 1 to n)
$4: FASTQ file containing the sample (integer, from 1 to n, matching to [FASTQ list file])

OPTIONS:
Option files contains two space-delimited columns: the optional argument is given in the first one and its value in the second one. Optional arguments are:
  - qual_threshold: Phred score for quality threshold averaged on 4 bp
  - min_sample_coverage: Minimum coverage per sample to keep a sequence
  - min_total_coverage: Minimum coverage over all sample to keep a sequence
  - min_median_coverage: COVmed, i.e. minimum median coverage of alleles
  - min_singleton_allele_coverage: COVsin, i.e. minimum coverage of a singleton allele
  - start_from_us: re-start analysis from existing .us files (1) or from FASTQ (0)
  - lane_lengths_given: number of reads given in [FASTQ list file] (1) or not (0)
  - min_ureads: minimum number of unique sequences required to keep a sample
  - save_library: save all intermediate files (1) or not (0)
  - interspecific_check: perform interspecific check (1) or not (0)
  - prior_sex_rate: prior sex rate in the population, in % (from 1 to 99)
  - max_nb_snps: maximum number of single nucleotide variation per sequence
  - max_nb_indels: maximum length of insertion-deletion per sequence
  

 

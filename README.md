# Allele-Express
Tools for Gene Expression Quantification in Polyploid Species
### 1、Dependencies

- ccs（ccs 6.4.0 (commit v6.4.0)）
- skera（skera 1.2.0）
- pbmm2  （pbmm2 1.13.1）
- samtools  （samtools 1.20  Using htslib 1.20）
- bedtools （bedtools v2.31.1）
- lima （lima 2.9.0）
- isoseq3  （isoseq 4.0.0 (commit v4.0.0)）
- isoquant  （IsoQuant 3.4.1）
- python 3.10.13 (pandas、matplotlib、numpy、seaborn、subprocess、configparser、argparse、string、collections、mpl_toolkits.axes_grid1.inset_locator、pathlib、random、pysam)
- R 4.4.0 (parallel、GenomicFeatures、tidyverse、data.table)
- mmseqs
- minimap2
- Allelefinder (For detailed usage instructions, please refer to https://github.com/sc-zhang/AlleleFinder)

### 2、 Installation and run

```shell
export PATH=/xx/scripts:$PATH
python run.py --config config.ini -o output
(After configuring config.ini, execute the command. Based on the number of threads set in config.ini, it will run in multiple threads.
Each step can independently adjust the number of threads.)
```

### 3、config.ini

- module_availability: true/false  Indicates whether this step needs to be executed. By default, subreads_to_hifi is false, AlleleFinder is false, and all others are true (please use absolute paths for the input data locations)。
  
  subreads_to_hifi: true/false (Control [0-ccs] module)
  
  read_segmentation: true/false (Control [1-read_segmentation] module)

  split_QC: true/false (Control [2-split_QC] module)

  demultiplex: true/false (Control [3-demultiplex] module)

  refine: true/false (Control [4-refine] module)

  ref_transcript_quant: true/false (Control [5-express] module)

  AlleleFinder: true/false (Control [6-AlleleFinder] module)

  Please select the appropriate step based on the data preprocessing situation and go ahead with the analysis.
  
- ref_data: Please input the reference CDS (in FA format), genome (in FA format), annotation (in standard GTF format, it is recommended to use AGAT (Another GTF/GFF Analysis Toolkit) software for format correction before running), and transcript (in FA format), respectively.

- [0-ccs]: Generate HIFI reads from subreads. This process is usually not employed. It depends on the type of sequencing data obtained.

  > 1. input_file: The input subreads file is in BAM format. For example: input_file = "subreads.bam"
  > 2. user_parameters： Parameters that the user needs to add themselves, such as: user_parameters = "--minZScore -5"
  > 3. The other parameters are identical to those of the CCS software.


- [1-read_segmentation]： Split the MAS connector

  > 1. primer_file： input the connection sequence of MAS. For example: primer_file ="/public/home/xx/mas8_primers.fasta"
  > 2. input_file： input the BAM file of "hifi reads", for example: input_file ="/public/home/xx/xx.bam"

- [2-split_QC]： Evaluate the results of the previous step of splitting

  > 1. cds_sample_number： The number of CDS samples. Example: cds_sample_number = 50000

- [3-demultiplex]： Remove the primer sequence

  > 1. primer_file2：Enter a primer sequence of 5'3'. Example: primer_file2 ="/public/home/xx/xx.primer.fasta"

- [4-refine]： Polymers with removed polyA tails and artificial primers

  > 1. The parameters are consistent with isoseq3 refine. The primer sequence to be entered here is the same as that entered in [3-demultiplex].

- [5-express]：allele-express section

  > 1. r_script_path = /public/home/xx/all_script
     It is the path of the R script. Because it has been called multiple times, it helps the script find the called R script
     
  > 2. sampling_proportions
     Proportional sampling is used to evaluate the comparison results. Example: "0.05 0.1 0.2 0.4 0.6 0.8 1.0"
     
  > 3. flnc_fa_file
     By default, it is empty. For details, please refer to the following precautions
     
  > 4. mmseqs2_env_path
    The path of mmseqs2, for example: mmseqs2_env_path = /xx/envs/mmseqs2
     
  > 5. delimiter
     Annotate the delimiters of polyploid homologous chromosomes in the first column of the GTF file. For example, there are four sets of homologous chromosomes in the gtf, and the first column is respectively:
       chr01_1、
       chr01_2、
       chr03_3、
       chr04_4、
       Then delimiter= _
     
- [6-AlleleFinder] is used to identify the alleles of polyploids. This step can run independently.

  > 1. polyploid_gff3 and polyploid_cds are annotations and cds data for which alleles need to be identified.

  > 2. ref_genome, ref_cds, and ref_gff3 are reference data for 8-AlleleFinder (usually a set of mono data).

  > 3. NUM_ALLELE： Ploidy of species

### 4、Important Notes 

Be careful about the formatting in the config.ini file. Some sections require double quotes, while others do not. Follow the format of the referenced config.ini file for adding. 

- [5-express]

flnc_fa_file: The flnc_fa_file is initially empty. If quantitative analysis is directly performed from this step, then it should be used as the input file. (If it is empty, the program will automatically search for the output file from [4-refine] as the input for this parameter) 

This step requires a very large amount of memory. If there is a memory error, it is recommended to reduce the number of threads in this step. 

- [6-AlleleFinder]

For detailed usage instructions, please refer to https://github.com/sc-zhang/AlleleFinder.

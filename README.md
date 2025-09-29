# Allele-Express
A novel algorithm specifically designed for allele-specific expression (ASE) analysis in haplotype-resolved genome assemblies.
## Table of contents
- [Prerequisites & Environment Setup](###prerequisites--environment-setup)
- [Getting Started](###getting-started)
- [Configuration Guide](###configuration-guide)
- [Usage Guidelines & Caveats](###usage-guidelines--caveats)
- [Output Files Documentation](###output-files-documentation)
## 1、Prerequisites & Environment Setup

- ccs（ccs 6.4.0 (commit v6.4.0)）
- skera（skera 1.2.0）
- pbmm2  （pbmm2 1.13.1）
- samtools  （samtools 1.20  Using htslib 1.20）
- bedtools （bedtools v2.31.1）
- lima （lima 2.9.0）
- isoseq3  （isoseq 4.0.0 (commit v4.0.0)）
- isoquant  （IsoQuant 3.4.1）**(optional)**
- python 3.10.13 （pandas、matplotlib、numpy、seaborn、subprocess、configparser、argparse、string、collections、mpl_toolkits.axes_grid1.inset_locator、pathlib、random、pysam）
- R 4.4.0 （parallel、GenomicFeatures、tidyverse、data.table）
- mmseqs （MMseqs2 Version: 15.6f452）
- minimap2 （minimap2 2.27-r1193）
- Allelefinder (For detailed usage instructions, please refer to https://github.com/sc-zhang/AlleleFinder)

## 2、Getting Started

```shell
#Resolve dependencies
export PATH =...
export PATH =/xx/scripts:$PATH
python run.py --config config.ini -o output
```
## 3、Configuration Guide  

Customize the application behavior by modifying the following parameters in config.ini:  

<ol>
<li> module_availability: true/false </li>  
  
This parameter indicates whether the step needs execution.  

By default, subreads_to_hifi and AlleleFinder are set to false, while all other steps are enabled.  
  
  > 1. subreads_to_hifi: true/false (Control [0-ccs] module)   

  > 2. read_segmentation: true/false (Control [1-read_segmentation] module)

  > 3. split_QC: true/false (Control [2-split_QC] module)   

  > 4. demultiplex: true/false (Control [3-demultiplex] module)   

  > 5. refine: true/false (Control [4-refine] module)   

  > 6. ref_transcript_quant: true/false (Control [5-express] module)   

  > 7. AlleleFinder: true/false (Control [6-AlleleFinder] module)


*Note: Absolute paths must be specified for input data locations. Please select the appropriate step based on the data preprocessing situation and go ahead with the analysis.*
  
<li>ref_data</li> 
The following files are required in specified formats:

  > 1. ref_cds: Format: FASTA (.fa or .fasta)
  > 2. ref_genome: Format: FASTA (.fa or .fasta)
  > 3. ref_gtf: Format: Standard GTF. *Recommendation: Pre-process the GTF file using AGAT (Another GTF/GFF Analysis Toolkit) to ensure format consistency.*
  > 5. ref_transcript: Format: FASTA (.fa or .fasta) *Selection Rule: For genes with multiple isoforms, the longest transcript per gene is used as the representative sequence to simplify analysis.*

<li>[0-ccs]</li> 
Generate HIFI reads from subreads. This process is usually not employed. It depends on the type of sequencing data obtained.

  > 1. input_file: The input subreads file is in BAM format. For example: input_file = "subreads.bam"
  > 2. user_parameters： Parameters that the user needs to add themselves, such as: user_parameters = "--minZScore -5"
  > 3. The other parameters are identical to those of the CCS software.


<li>[1-read_segmentation]</li> 
Split the MAS connector

  > 1. primer_file： input the connection sequence of MAS. For example: primer_file ="/public/home/xx/mas8_primers.fasta"
  > 2. input_file： input the BAM file of "hifi reads", for example: input_file ="/public/home/xx/xx.bam"

<li>[2-split_QC]</li> 
Evaluate the results of the previous step of splitting

  > 1. cds_sample_number： The number of CDS samples. Example: cds_sample_number = 50000

<li>[3-demultiplex]</li> 
Remove the primer sequence

  > 1. primer_file2：Enter a primer sequence of 5'3'. Example: primer_file2 ="/public/home/xx/xx.primer.fasta"

<li>[4-refine]</li> 
Polymers with removed polyA tails and artificial primers

  > 1. The parameters are consistent with isoseq3 refine. The primer sequence to be entered here is the same as that entered in [3-demultiplex].

<li>[5-express]</li>
allele-express section

  > 1. r_script_path = /public/home/xx/all_script
     It is the path of the R script. Because it has been called multiple times, it helps the script find the called R script
     
  > 2. sampling_proportions
     Proportional sampling is used to evaluate the comparison results. Example: "0.05 0.1 0.2 0.4 0.6 0.8 1.0"
     
  > 3. flnc_fa_file
     By default, it is empty. For details, please refer to the following precautions
     
  > 4. mmseqs2_env_path
    The path of mmseqs2, for example: mmseqs2_env_path = /xx/envs/mmseqs2
     
  > 5. delimiter
     Set the delimiter (e.g., '_') according to the chromosome naming pattern (e.g., 'chr01_1') in column 1 of the GTF file.:
     
<li>[6-AlleleFinder]</li>
It is used to identify the alleles of polyploids. This step can run independently.

  > 1. polyploid_gff3 and polyploid_cds are annotations and cds data for which alleles need to be identified.

  > 2. ref_genome, ref_cds, and ref_gff3 are reference data for 8-AlleleFinder (usually a set of mono data).

  > 3. NUM_ALLELE： Ploidy of species
  </ol>
  
## 4、Usage Guidelines & Caveats
Before deploying, please review these critical constraints and best practices:
Be careful about the formatting in the config.ini file. Some sections require double quotes, while others do not. Follow the format of the referenced config.ini file for adding. 
<ol>
<li>[5-express]</li>

- flnc_fa_file: The flnc_fa_file is initially empty. If quantitative analysis is directly performed from this step, then it should be used as the input file. (If it is empty, the program will automatically search for the output file from [4-refine] as the input for this parameter) 

- threads: This step requires a very large amount of memory. If there is a memory error, it is recommended to reduce the number of threads in this step. (The recommended number of threads should not exceed 10.)

<li>[6-AlleleFinder]</li>

- For detailed usage instructions, please refer to https://github.com/sc-zhang/AlleleFinder.
</ol>

## 5、Output Files Documentation
This pipeline generates the following result files for transcript-level and gene-level quantification:  

<ol>
<li>File List</li>

- trans2gene.txt
- trans2trans.txt
- transcript_counts.txt
- gene_counts.txt

<li>trans2gene.txt</li>
Purpose: Maps transcripts to their corresponding genes, serving as the basis for gene-level quantification.  

Format:

`TXNAME    GENEID    Length`

<li>trans2trans.txt</li>
Purpose: Identifies transcript groups with identical sequences and designates a representative transcript for quantification to avoid redundancy.  

Format:

`representative_transcript    transcript_group`

Example:

`POJ2078104349.t1    POJ2078104349.t1,POJ2078104345.t1`

Note:
- Transcripts not selected as representatives (e.g., POJ2078104345.t1) will not appear directly in transcript_counts.txt.
  Their expression values are derived from the representative transcript (e.g., POJ2078104349.t1).
  
- The gene-level quantification (gene_counts.txt) follows the same principle: each gene's expression value is the sum of counts from all its representative transcripts (mapped via trans2gene.txt).
  
- To verify non-representative transcripts, check their mappings in this file.

<li>transcript_counts.txt</li>
Purpose: Records raw expression counts (Counts) for representative transcripts.  

Format:

`TranscriptID    Counts`

<li>gene_counts.txt</li>
Purpose: Aggregates gene-level expression data, including raw counts and normalized values (TPM).  

Format:

`GeneID    Counts    TPM`

</ol>

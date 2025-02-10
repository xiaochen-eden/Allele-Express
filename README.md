# Allele-Express
Tools for Gene Expression Quantification in Polyploid Species
### 1、依赖

- ccs（ccs 6.4.0 (commit v6.4.0)）
- skera（skera 1.2.0）
- pbmm2  （pbmm2 1.13.1）
- samtools  （samtools 1.20  Using htslib 1.20）
- bedtools （bedtools v2.31.1）
- lima （lima 2.9.0）
- isoseq3  （isoseq 4.0.0 (commit v4.0.0)）
- isoquant  （IsoQuant 3.4.1）
- python 3.10.13(pandas、matplotlib、numpy、seaborn、subprocess、configparser、argparse、string、collections、mpl_toolkits.axes_grid1.inset_locator、pathlib、random、pysam)
- R 4.4.0(parallel、GenomicFeatures、tidyverse、data.table)
- mmseqs
- minimap2
- Allelefinder及其依赖(具体使用教程参照https://github.com/sc-zhang/AlleleFinder)

### 2、安装与运行

```shell
export PATH=/xx/scripts:$PATH
python run.py --config config.ini -o output  （配置完config.ini之后提交命令，根据contig.ini中设置的线程数多线程运行,可以每个步骤单独调整线程数）
```

### 3、config.ini

- module_availability：true/false  表示是否需要运行该步骤，默认情况下，subreads_to_hifi为false,其他均为true（输入数据的位置请均使用绝对路径）
- ref_data 请分别输入参考cds（fa格式）、基因组（fa格式）、注释（gtf格式）、转录本（fa格式）

- [0-ccs]: 将subreads生成hifi reads,一般不使用，根据下机数据类型判断

  > 1. input_file : 输入的subreads文件，bam格式。比如：input_file = "subreads.bam"
  > 2. user_parameters：用户自己需要添加的参数，比如：user_parameters = "--minZScore -5"
  > 3. 其他参数与ccs软件的参数一致。


- [1-read_segmentation]：拆分MAS接头

> 1. primer_file：输入MAS的接头序列  比如：primer_file ="/public/home/xx/mas8_primers.fasta"
> 2. input_file：输入hifi reads的bam文件  比如：input_file ="/public/home/xx/xx.bam"

- [2-split_QC]：对上一步拆分结果进行评估

> cds_sample_number：抽样的cds个数。例子：cds_sample_number = 50000

- [3-demultiplex]：去除引物序列

> primer_file2：输入5'3'的引物序列。例子：primer_file2 ="/public/home/xx/xx.primer.fasta"

- [4-refine]：去除polyA尾和人工引物的多聚体

> 参数与isoseq3 refine一致，此处需输入的引物序列与[3-demultiplex]输入的引物序列一致。

- [5-express]：allele-express部分

  1. r_script_path = /public/home/xx/all_script
     是R脚本的路径，因为执行了多次调用，帮助脚本找到调用的R脚本
     
  2. sampling_proportions
     按比例抽样用来评估比对结果。例子："0.05 0.1 0.2 0.4 0.6 0.8 1.0"
     
  3. flnc_fa_file
     默认为空，详见下面注意事项
     
  4. mmseqs2_env_path
     mmseqs2的路径，例子：mmseqs2_env_path = /xx/envs/mmseqs2
     
  5. delimiter
     注释gtf文件中第一列的多倍体同源染色体的分割符，例子，gtf中有四套同源染色体，第一列分别为chr01_1 chr01_2 chr03_3 chr04_4，则delimiter= _
     
- [6-AlleleFinder]

> assembly_gff3、assembly_cds、assembly_genome 为需要鉴定等位的基因组数据。

> 要注意[ref_data]模块中，ref_cds与ref_gff3为8-AlleleFinder的参考数据（一般为一套mono数据）。

> NUM_ALLELE：倍性

### 4、注意事项

config.ini中注意格式问题，有的需要加双引号，有的不需要加，按照参考的config.ini格式添加。

[5-GQ-mapping]

flnc_fq_file这里是如果前面的数据已经通过别的软件处理好，整个流程从5-mapping开始，则flnc_fq_file作为输入文件，从此处开始。（前提是[4-refine]中并没有处理好的reads数据）

[7-express]

flnc_fa_file这里是如果前面的数据已经通过别的软件处理好,整个流程从7-express开始,则flnc_fa_file作为输入文件，从此处开始。(前提是[4-refine]中并没有处理好的reads数据)
该步骤需要内存非常大，如果出现内存报错，建议降低该步骤的线程数。

[8-AlleleFinder]

具体使用教程参照https://github.com/sc-zhang/AlleleFinder。

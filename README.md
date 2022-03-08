# gatk pipeline based on python
----------------------------------------------
gatk calling variants include two mode, the first one is single calling (produce gvcf for each sample), the second one is population calling (directly produce vcf include all  sample).<br/>

It also supports genome splitting by setting -ss or -sn, if you split the genome, it will calling variant on each split region.<br/>

specify the input directory by setting the -i parameter, It will automatically look for files ending with suffix_name (by setting the -S options)in this directory.


## Usage:
```
 usage: gatk.py [-h] [-v] -i  -o  [-m] [-sn  | -ss ] [--split-mode] [-hcjo] [-cgjo] [-ggjo] [-mvjo] [--gvcf] [-s  [...]] [-p] [-pbi] [-t] -r  [-n] [-l]

The purpose of this program is for calling snp used gatk!
--------------------------------------------------------------------

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -i , --indir          the directory where bam located
  -o , --outdir         the directory save the snp and indel calling result
  -m , --mode           the calling mode, pc: population calling, sc: single calling [sc]
  -sn , --split-num     the number of split genome
  -ss , --split-size    the size for split genome [unit M]
  --split-mode          the mode of split genome [imprecise]
  -hcjo , --haplotypecaller-java-options
                        the java options for HaplotypeCaller, such as -Xmx100G limit the memory usage
  -cgjo , --combinegvcf-java-options
                        the java options for CombineGVCFs
  -ggjo , --genotypegvcf-java-options
                        the java options for GenotypeGVCFs
  -mvjo , --mergevcf-java-options
                        the java options for MergeVcfs
  --gvcf                output gvcf file , default output is vcf file [False]
  -s  [ ...], --suffix-name  [ ...]
                        the suffix name of input file [default sort.fixmate.rmdup.bam]
  -p , --process        the porcess number for deal with multiple samples [5]
  -pbi , --process-build-index
                        the process of build bam index [5]
  -t , --threads        the threads for index and sort bam file [5]
  -r , --ref            the reference genome
  -n , --prefix-name    the prefix name of output vcf file [All]
  -l , --log            the log file for record running information [gatk.log]

For command line option of each command, type: gatk.py COMMAND -h



# Sigma: Strain-level Identification of Genomes from Metagenomic Analysis for Biosurveillance

### Cite this article
Tae-Hyuk Ahn, Juanjuan Chai, Chongle Pan, Sigma: Strain-level inference of genomes from metagenomic analysis for biosurveillance, Bioinformatics, Volume 31, Issue 2, January 2015, Pages 170–177, [https://doi.org/10.1093/bioinformatics/btu641](https://doi.org/10.1093/bioinformatics/btu641)

### Summary
Identification of microbes and accurate estimation of their relative abundance are crucial subjects in metagenomic analysis. Existing taxonomic classification methods are unsuitable for metagenomic biosurveillance due to the following three factors. First, reference genomes including complete genomes of pathogenic bacteria are required. Second, strain-level genome identification from the same species is essential. Third, statistical confidence evaluation of the metagenomic analysis is recommended. We developed Sigma (Strain-level Inference of Genomes from Metagenomic Analysis), a novel sequence similarity-based approach for strain-level identification of genomes from metagenomic analysis for biosurveillance. Nucleotide-level alignment analysis using maximum likelihood estimation empowers Sigma to estimate the relative abundances and likelihood ratios of each genome accurately given a list of reference genomes. The hybrid parallel architecture of the SIGMA allows its computation to be scalable from desktops to supercomputers for different sizes of metagenomic reads and reference datasets.
![](https://github.com/BioHPC/Sigma/blob/main/figures/overview_figure.gif){:style="display:block; margin-left:auto; margin-right:auto"}

### Manual
Please check the SIGMA User Manual for details.

### Prerequisites
* Bowtie2 [https://github.com/BenLangmead/bowtie2](https://github.com/BenLangmead/bowtie2)
* Samtools [http://www.htslib.org/](http://www.htslib.org/)

### Install and Run Sigma
Binaries were built and tested in the Intel architectures (x86_64) running Linux/Unix environment. Sigma package was developed using C++. C++ source codes were compiled and linked by g++ (gcc version 4.7.2).

Please download the last version of Sigma and extract it from the command line

$ tar xzvf sigma-version.tar.gz
After you extract the tar ball, please check usage with –h option. For example,

$ cd Sigma/bin
$ ./sigma -h
Please read the user manual to run the Sigma package.

### Sigma Algorithm
The foundation of all Sigma’s statistical inferences is a simple probabilistic model for sampling reads with mismatches from genomes at unknown abundances. The Sigma probabilistic model permits maximum likelihood estimation (MLE) of the relative abundances of all genomes in terms of the percentages of reads sampled from these genomes. Importantly, we proved the convexity of the objective function of the MLE such that Sigma can guarantee to attain a global optimum solution using interior point non-linear programming (NLP) provided in the Ipopt library.

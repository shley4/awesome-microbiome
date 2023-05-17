# awesome_meta-omics

[![Awesome](https://awesome.re/badge.svg)](https://awesome.re)

List of software/packages for meta-omic data analysis, including assembly (short reads, long reads, and hybrid), genome recovery, annotation, taxonomy, abandance, and so on. [Contributions welcome](https://github.com/goekelab/awesome-nanopore#contributing)...
## sequencing data quality control
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - [?] - A quality control tool for high throughput sequence data.
- [fastp](https://github.com/OpenGene/fastp) - [C++/C] - [fastp: an ultra-fast all-in-one FASTQ preprocessor](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234)
- [Trimmomatic](https://github.com/usadellab/Trimmomatic) - [Java] - [Trimmomatic: a flexible trimmer for Illumina sequence data](https://academic.oup.com/bioinformatics/article/30/15/2114/2390096)
- [Cutadapt](https://github.com/marcelm/cutadapt) - [Python] - [Cutadapt Removes Adapter Sequences From High-Throughput Sequencing Reads](https://doi.org/10.14806/ej.17.1.200)
- [FastUniq](https://sourceforge.net/projects/fastuniq/files/) - [C++?] - [FastUniq: A Fast De Novo Duplicates Removal Tool for Paired Short Reads](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0052249)

## Amplicon
- [EasyAmplicon](https://github.com/YongxinLiu/EasyAmplicon) - [EasyAmplicon: An easy-to-use, open-source, reproducible, and community-based pipeline for amplicon data analysis in microbiome research](https://onlinelibrary.wiley.com/doi/full/10.1002/imt2.83). Beautiful figures.
- [dada2](https://benjjneb.github.io/dada2/index.html) - [R] [DADA2: High resolution sample inference from Illumina amplicon data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/)
- [QIIME2](https://qiime2.org/) - [?] - [Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2](https://doi.org/10.1038/s41587-019-0209-9)

## assembly
### assembly - short reads
- [metaSPAdes](https://github.com/ablab/spades) - [C++] - [metaSPAdes: a new versatile metagenomic assembler](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5411777/)
- [megahit](https://github.com/voutcn/megahit) - [C++] - Ultra-fast and memory-efficient short read assembler, but also works well on single genome and single-cell assembly. - [MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph](https://academic.oup.com/bioinformatics/article/31/10/1674/177884)

### assembly - long reads
- [metaFlye](https://github.com/fenderglass/Flye) - [C/C++] - [metaFlye: scalable long-read metagenome assembly using repeat graphs](metaFlye: scalable long-read metagenome assembly using repeat graphs)
- [Canu](https://github.com/marbl/canu) - [C++/C] - [Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation](https://genome.cshlp.org/content/27/5/722)
- [FALCON](https://github.com/PacificBiosciences/FALCON) - [Pytohn]
- [miniasm](https://github.com/lh3/miniasm) - [TeX/C] - [Minimap and miniasm: fast mapping and de novo assembly for noisy long sequences](https://academic.oup.com/bioinformatics/article/32/14/2103/1742895)
- [wtdbg2](https://github.com/ruanjue/wtdbg2) - [C] - [Fast and accurate long-read assembly with wtdbg2](https://www.nature.com/articles/s41592-019-0669-3)

### assembly - hybrid
- [metaSPAdes](https://github.com/ablab/spades) - [C++] - metagenomic assembly
- [Aviary](https://github.com/rhysnewell/aviary) - [Python]
- [OPERA-MS](https://github.com/CSB5/OPERA-MS) - [C++/Python] - [Hybrid metagenomic assembly enables high-resolution analysis of resistance determinants and mobile elements in human microbiomes](https://www.nature.com/articles/s41587-019-0191-2)


## genomes recovery
### tools and workflows
- [MetaBAT2](https://bitbucket.org/berkeleylab/metabat) - [C++] - [MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6662567/)
- [MaxBin2](https://sourceforge.net/projects/maxbin2/) - [C++] - [MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets](https://academic.oup.com/bioinformatics/article/32/4/605/1744462)
- [CONCOCT](https://github.com/BinPro/CONCOCT) - [Python] - [Binning metagenomic contigs by coverage and composition](https://www.nature.com/articles/nmeth.3103)
- [VAMB](https://github.com/RasmussenLab/vamb) - [Python] - [Improved metagenome binning and assembly using deep variational autoencoders](https://www.nature.com/articles/s41587-020-00777-4)
- [UniteM](https://github.com/donovan-h-parks/UniteM) - [Python] - a software toolkit implementing different ensemble binning strategies for producing a non-redundant set of bins from the output of multiple binning methods.

## read mapping
### short-read mapping
- [Bowtie2](https://github.com/BenLangmead/bowtie2) - [C++] - [Fast gapped-read alignment with Bowtie 2](https://www.nature.com/articles/nmeth.1923)
- [bwa](https://github.com/lh3/bwa) - [C] - [Fast and accurate short read alignment with Burrows-Wheeler transform](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2705234/)

### long-read mapping
- [minimap2](https://github.com/lh3/minimap2) - [C] - [Minimap2: pairwise alignment for nucleotide sequences](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778)

## genome abundance estimation
- [coverM](https://github.com/wwood/CoverM) - [Rust] -  fast DNA read coverage and relative abundance calculator focused on metagenomics applications.

## genome assessment
- [CheckM](https://github.com/Ecogenomics/CheckM) - [Python] - [CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes](https://pubmed.ncbi.nlm.nih.gov/25977477/)
- [CheckM2](https://github.com/chklovski/CheckM2) - [Scilab] - [bioRxiv - CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning](https://www.biorxiv.org/content/10.1101/2022.07.11.499243v1)
- [BUSCO](https://gitlab.com/ezlab/busco) - [arXiv- BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes](http://arxiv.org/abs/2106.11799)
- [GTDBTk](https://github.com/Ecogenomics/GTDBTk) - [Python] - [ GTDB-Tk v2: memory friendly classification with the Genome Taxonomy Database](https://academic.oup.com/bioinformatics/article/38/23/5315/6758240?login=false&utm_source=advanceaccess&utm_campaign=bioinformatics&utm_medium=email)
- [RefineM](https://github.com/donovan-h-parks/RefineM) - [Python] - A toolbox for improving metagenome-assembled genomes.
- [MAGpurify](https://github.com/snayfach/MAGpurify) - [Python] - [New insights from uncultivated genomes of the global human gut microbiome](https://www.nature.com/articles/s41586-019-1058-x)


## gene prediction
- [Prodigal](https://github.com/hyattpd/Prodigal) - [C] - [Prodigal: prokaryotic gene recognition and translation initiation site identification](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119)
- [Prokka](https://github.com/tseemann/prokka) - [Perl] - [Prokka: rapid prokaryotic genome annotation](http://www.ncbi.nlm.nih.gov/pubmed/24642063)

## gene annotation
### tools
- [eggNOG-mapper v2](https://github.com/eggnogdb/eggnog-mapper) - [Python] - [eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale](https://doi.org/10.1093/molbev/msab293)
- [KofamScan](https://www.genome.jp/ftp/tools/kofam_scan/) - [NA] - [KofamKOALA: KEGG ortholog assignment based on profile HMM and adaptive score threshold](https://doi.org/10.1093/bioinformatics/btz859)

### databases
- [eggNOG 6.0](http://eggnog6.embl.de/) - [eggNOG 6.0: Enabling comparative genomics across 12,535 organisms.](https://doi.org/10.1093/nar/gkac1022)
- [KEGG]

## phylogenetics analysis
### phylogeny construction
- [tree]
- []

### phylogeny visualization
- [iTol]
- [ARB]

## comparative genomics
- [fastANI]
- [AAI]
- [synteny]
- [snp]
- [indel]
- [sv]

## workflows
- [Aviary](https://github.com/rhysnewell/aviary) - [Python] - a snakemake pipeline for assembly, binning, annotation, and strain diversity analysis.
- [MetaWRAP](https://github.com/bxlab/metaWRAP) - [Shell/Python] - [MetaWRAP—a flexible pipeline for genome-resolved metagenomic data analysis](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0541-1)


## Contributing
We welcome contributions and suggestions! Please follow the steps below to contribute:
1. Fork this repository
2. Make a change to README.md in this format: `[RESOURCE](LINK)` - [language(s)] - DESCRIPTION
3. Submit a pull request.

## Citation
Wan, Y.K., Hendra, C., Pratanwanich, P.N. & Göke, J. Beyond sequencing: machine learning algorithms extract biology hidden in Nanopore signal data. *Trends in Genetics* (2021). https://doi.org/10.1016/j.tig.2021.09.001.

## What is an awesome list?
According to [the official awesome Github repository](https://github.com/sindresorhus/awesome), an awesome list on GitHub is "a curation of actual awesome stuff", so an awesome list only includes items that has been reseached on by a contributor who would personally recommend the items. To learn more, please read [the official awesome manifesto](https://github.com/sindresorhus/awesome/blob/main/awesome.md).

## Contact
This repository was created by [Jie Li](https://github.com/jlli6t).

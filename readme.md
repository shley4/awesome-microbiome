# magic_microbiome

[![Awesome](https://awesome.re/badge.svg)](https://awesome.re)

List of tools for microbiome data.

## sequencing data quality control
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - [?] - A quality control tool for high throughput sequence data.
- [fastp](https://github.com/OpenGene/fastp) - [C++/C] - [fastp: an ultra-fast all-in-one FASTQ preprocessor](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234)
- [Trimmomatic](https://github.com/usadellab/Trimmomatic) - [Java] - [Trimmomatic: a flexible trimmer for Illumina sequence data](https://academic.oup.com/bioinformatics/article/30/15/2114/2390096)
- [Cutadapt](https://github.com/marcelm/cutadapt) - [Python] - [Cutadapt Removes Adapter Sequences From High-Throughput Sequencing Reads](https://doi.org/10.14806/ej.17.1.200)
- [FastUniq](https://sourceforge.net/projects/fastuniq/files/) - [C++] - [FastUniq: A Fast De Novo Duplicates Removal Tool for Paired Short Reads](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0052249)
- [PoreChop][https://github.com/rrwick/Porechop] - [C++] - adapter trimmer for Oxford Nanopore reads - Unsupported since Oct 2018

## polish long reads or assemblies
- [FMLRC2](https://github.com/HudsonAlpha/fmlrc2) - [Rust] - [FMLRC2 reference: Mak, Q. C., _et al_. (2023). Polishing de novo nanopore assemblies of bacteria and eukaryotes with FMLRC2. _Molecular Biology and Evolution_, 40(3), msad048](https://doi.org/10.1093/molbev/msad048)
- [pilon](https://github.com/broadinstitute/pilon) - [Scala] - [Bruce J. Walker, _et al_. (2014) Pilon: An Integrated Tool for Comprehensive Microbial Variant Detection and Genome Assembly Improvement. _PLoS ONE_ 9(11): e112963](https://doi.org/10.1371/journal.pone.0112963)
- [NextPolish](https://github.com/Nextomics/NextPolish) - [NA] - [Hu, Jiang, _et al_. “NextPolish: a fast and efficient genome polishing tool for long read assembly.” Bioinformatics (Oxford, England) (2019)](https://doi.org/10.1093/bioinformatics/btz891)
- [Polypolish](https://github.com/rrwick/Polypolish) - [Rust] - [Bouras G, Judd LM, Edwards RA, Vreugde S, Stinear TP, Wick RR. How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies. _Microbial Genomics_. 2024](https://doi.org/10.1099/mgen.0.001254)
- [pypolca](https://github.com/gbouras13/pypolca) - [Python] - [Bouras G, Judd LM, Edwards RA, Vreugde S, Stinear TP, Wick RR. How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies. _Microbial Genomics_. 2024](https://doi.org/10.1099/mgen.0.001254)

## Amplicon
### tools
- [dada2](https://benjjneb.github.io/dada2/index.html) - [R] [DADA2: High resolution sample inference from Illumina amplicon data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/)
- [QIIME2](https://qiime2.org/) - [NA] - [Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2](https://doi.org/10.1038/s41587-019-0209-9)
[RDP Classifier](https://sourceforge.net/projects/rdp-classifier/) - [Updated RDP taxonomy and RDP Classifier for more accurate taxonomic classification](https://journals.asm.org/doi/10.1128/mra.01063-23)

### databases
- [EUKARYOME](https://eukaryome.org) - [EUKARYOME: the rRNA gene reference database for identification of all eukaryotes](https://doi.org/10.1093/database/baae043)
- [Silva](https://www.arb-silva.de) - [The SILVA ribosomal RNA gene database project: improved data processing and web-based tools](http://nar.oxfordjournals.org/content/41/D1/D590)
- [Greengenes2](http://ftp.microbio.me/greengenes_release/) - [Greengenes2 unifies microbial data in a single reference tree](https://doi.org/10.1038/s41587-023-01845-1)
- [MiDAS](https://midasfieldguide.org/guide) - [Dueholm, M.K.D., Andersen, K.S., Korntved, AK.C. _et al_. MiDAS 5: Global diversity of bacteria and archaea in anaerobic digesters. _Nat Commun_ 15, 5361 (2024)](https://doi.org/10.1038/s41467-024-49641-y)

## assembly
### assembly - short reads
- [metaSPAdes](https://github.com/ablab/spades) - [C++] - [metaSPAdes: a new versatile metagenomic assembler](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5411777/)
- [megahit](https://github.com/voutcn/megahit) - [C++] - [MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph](https://academic.oup.com/bioinformatics/article/31/10/1674/177884)
- [PenuiN](https://github.com/soedinglab/plass) - [C] - [Annika Jochheim, _et al_., _bioRxiv_, 2024. Strain-resolved de-novo metagenomic assembly of viral genomes and microbial 16S rRNAs](https://doi.org/10.1101/2024.03.29.587318)

### assembly - long reads
- [metaFlye](https://github.com/fenderglass/Flye) - [C/C++] - [metaFlye: scalable long-read metagenome assembly using repeat graphs](https://doi.org/10.1038/s41592-020-00971-x)
- [Canu](https://github.com/marbl/canu) - [C++/C] - [Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation](https://genome.cshlp.org/content/27/5/722)
- [FALCON](https://github.com/PacificBiosciences/FALCON) - [Pytohn]
- [miniasm](https://github.com/lh3/miniasm) - [TeX/C] - [Minimap and miniasm: fast mapping and de novo assembly for noisy long sequences](https://academic.oup.com/bioinformatics/article/32/14/2103/1742895)
- [wtdbg2](https://github.com/ruanjue/wtdbg2) - [C] - [Fast and accurate long-read assembly with wtdbg2](https://www.nature.com/articles/s41592-019-0669-3)
- [Trycycler](https://github.com/rrwick/Trycycler) - [Pyton] - [Wick RR, _et al_. Trycycler: consensus long-read assemblies for bacterial genomes. _Genome Biology_. 2021](https://doi.org/10.1186/s13059-021-02483-z)
- [NextDenovo](https://github.com/Nextomics/NextDenovo) - [C] - [Hu, J., Wang, Z., Sun, Z. _et al_. NextDenovo: an efficient error correction and accurate assembly tool for noisy long reads. _Genome Biol_ 25, 107 (2024)](https://doi.org/10.1186/s13059-024-03252-4)
- [hifiasm-meta](https://github.com/xfengnefx/hifiasm-meta/) - [C++] - [Feng, X., Li, H. Evaluating and improving the representation of bacterial contents in long-read metagenome assemblies. _Genome Biol_ 25, 92 (2024)](https://doi.org/10.1186/s13059-024-03234-6)

### assembly - hybrid
- [metaSPAdes](https://github.com/ablab/spades) - [C++]
- [Aviary](https://github.com/rhysnewell/aviary) - [Python]
- [OPERA-MS](https://github.com/CSB5/OPERA-MS) - [C++/Python] - [Hybrid metagenomic assembly enables high-resolution analysis of resistance determinants and mobile elements in human microbiomes](https://www.nature.com/articles/s41587-019-0191-2)
- [Unicycler](https://github.com/rrwick/Unicycler) - [C++ ] - [Wick RR, Judd LM, Gorrie CL, Holt KE. Unicycler: resolving bacterial genome assemblies from short and long sequencing reads. PLoS Comput Biol 2017](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005595)

## binning
### tools and workflows
- [MetaBAT2](https://bitbucket.org/berkeleylab/metabat) - [C++] - [MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6662567/)
- [MaxBin2](https://sourceforge.net/projects/maxbin2/) - [C++] - [MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets](https://academic.oup.com/bioinformatics/article/32/4/605/1744462)
- [CONCOCT](https://github.com/BinPro/CONCOCT) - [Python] - [Binning metagenomic contigs by coverage and composition](https://www.nature.com/articles/nmeth.3103)
- [VAMB](https://github.com/RasmussenLab/vamb) - [Python] - [Improved metagenome binning and assembly using deep variational autoencoders](https://www.nature.com/articles/s41587-020-00777-4)
- [UniteM](https://github.com/dparks1134/UniteM) - [Python]
- [SemiBin](https://github.com/BigDataBiology/SemiBin) - [Python] - [A deep siamese neural network improves metagenome-assembled genomes in microbiome datasets across different environments](https://doi.org/10.1038/s41467-022-29843-y)
- [SemiBin2](https://github.com/BigDataBiology/SemiBin) - [Python] - [SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing](https://doi.org/10.1093/bioinformatics/btad209) - [Benchmark](https://github.com/BigDataBiology/SemiBin2_benchmark)
- [Smeta](https://github.com/YuhaoZhangwow/SMeta) - [C++] - [Yuhao Zhang, _et al_., _bioRxiv_, 2024. SMeta, a binning tool using single-cell sequences to aid in reconstructing species from metagenome accurately](https://doi.org/10.1101/2024.08.25.609542)
- [Bin Chiken](https://github.com/AroneyS/binchicken) - [Python]

## Taxonomy
### for mag
- [GTDBTk](https://github.com/Ecogenomics/GTDBTk) - [Python] - [ GTDB-Tk v2: memory friendly classification with the Genome Taxonomy Database](https://doi.org/10.1093/bioinformatics/btac672)

### for read
- [MetaPhlAn](https://github.com/biobakery/MetaPhlAn) - [Python] - [Aitor Blanco-Miguez, _et al_. Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. _Nature Biotechnology (2023)](https://doi.org/10.1038/s41587-023-01688-w)
- [StrainPhlAn](https://github.com/biobakery/MetaPhlAn) - [Python] - [Duy Tin Truong, _et al_. Microbial strain-level population structure and genetic diversity from metagenomes. _Genome Research_ (2017)](https://dx.doi.org/10.1101/gr.216242.116)
- [kraken2](https://github.com/DerrickWood/kraken2) - [C++] - [Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. _Genome Biol_ 20, 257 (2019)](https://doi.org/10.1186/s13059-019-1891-0)
- [Metabuli](https://github.com/steineggerlab/Metabuli) - [C++] - [Kim, J., Steinegger, M. Metabuli: sensitive and specific metagenomic classification via joint analysis of amino acid and DNA. _Nat Methods_ 21, 971–973 (2024)](https://doi.org/10.1038/s41592-024-02273-y)

## alignment
### short-read to sequence
- [Bowtie2](https://github.com/BenLangmead/bowtie2) - [C++] - [Fast gapped-read alignment with Bowtie 2](https://www.nature.com/articles/nmeth.1923)
- [bwa](https://github.com/lh3/bwa) - [C] - [Fast and accurate short read alignment with Burrows-Wheeler transform](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2705234/)

### long-read to sequence
- [minimap2](https://github.com/lh3/minimap2) - [C] - [Minimap2: pairwise alignment for nucleotide sequences](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778)

### sequence to sequence
- [blast+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) - [NA] - [Christiam Camacho, _et al_. BLAST+: architecture and applications. _BMC Bioinformatics_. (2009)](https://doi.org/10.1186%2F1471-2105-10-421) - [download BLAST databases](https://ftp.ncbi.nlm.nih.gov/blast/db/)
- [LexicMap](https://github.com/shenwei356/LexicMap) - [Go] [Wei Shen, Zamin Iqbal. 2024. LexicMap: efficient sequence alignment against millions of prokaryotic genomes](https://doi.org/10.1101/2024.08.30.610459) - bioRxiv
- [foldmason](https://github.com/steineggerlab/foldmason) - [Cameron Gilchrist, _et al_., _bioRxiv_, 2024. Multiple Protein Structure Alignment at Scale with FoldMason](https://doi.org/10.1101/2024.08.01.606130)

## abundance estimation
- [coverM](https://github.com/wwood/CoverM) - [Rust] -  fast DNA read coverage and relative abundance calculator focused on metagenomics applications.
- [Fairy](https://github.com/bluenote-1577/fairy) - [Rust] - [Shaw, J., Yu, Y. Fairy: fast approximate coverage for multi-sample metagenomic binning. _Microbiome_ 12, 151 (2024)](https://doi.org/10.1186/s40168-024-01861-6)
- [SingleM Microbial Fraction(SMF)](https://github.com/EisenRa/2024_soil_dark_matter_reply) - [Raphael Eisenhofer, Antton Alberdi, Ben J Woodcroft, Quantifying microbial DNA in metagenomes improves microbial trait estimation, *ISME Communications*, 2024](https://doi.org/10.1093/ismeco/ycae111)

## mag assessment
- [BUSCO](https://gitlab.com/ezlab/busco) - [arXiv- BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes](http://arxiv.org/abs/2106.11799)
- [CheckM](https://github.com/Ecogenomics/CheckM) - [Python] - [CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes](https://pubmed.ncbi.nlm.nih.gov/25977477/)
- [CheckM2](https://github.com/chklovski/CheckM2) - [Scilab] - [bioRxiv - CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning](https://www.biorxiv.org/content/10.1101/2022.07.11.499243v1)
- [RefineM](https://github.com/donovan-h-parks/RefineM) - [Python] - A toolbox for improving metagenome-assembled genomes.
- [MAGpurify](https://github.com/snayfach/MAGpurify) - [Python] - [New insights from uncultivated genomes of the global human gut microbiome](https://www.nature.com/articles/s41586-019-1058-x)
[DFAST_QC](https://github.com/nigyta/dfast_qc) - [Python] - [DFAST_QC: Quality Assessment and Taxonomic Identification Tool for Prokaryotic Genomes](https://doi.org/10.1101/2024.07.22.604526) - bioRxiv

## gene prediction
- [Prodigal](https://github.com/hyattpd/Prodigal) - [C] - [Prodigal: prokaryotic gene recognition and translation initiation site identification](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119)
- [Prokka](https://github.com/tseemann/prokka) - [Perl] - [Prokka: rapid prokaryotic genome annotation](http://www.ncbi.nlm.nih.gov/pubmed/24642063)

## structure
- [Chai-1](https://github.com/chaidiscovery/chai-lab) - [Python] - [Introducing Chai-1: Decoding the molecular interactions of life](https://www.chaidiscovery.com/blog/introducing-chai-1)

## annotation
### tools
- [eggNOG-mapper v2](https://github.com/eggnogdb/eggnog-mapper) - [Python] - [eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale](https://doi.org/10.1093/molbev/msab293)
- [KofamScan](https://www.genome.jp/ftp/tools/kofam_scan/) - [NA] - [KofamKOALA: KEGG ortholog assignment based on profile HMM and adaptive score threshold](https://doi.org/10.1093/bioinformatics/btz859)
- [Macrel](https://github.com/BigDataBiology/macrel) - [Python] - [Macrel: antimicrobial peptide screening in genomes and metagenomes](https://doi.org/10.7717/peerj.10555)
- [RGI](https://card.mcmaster.ca/analyze) - [Python] - [CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database](https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/36263822/)
- [MetaCerberus](https://github.com/raw-lab/metacerberus) - [Python] - [Jose L Figueroa III, _et al_.  MetaCerberus: distributed highly parallelized HMM-based processing for robust functional annotation across the tree of life. _Bioinformatics_ 40, 2024, btae119](https://doi.org/10.1093/bioinformatics/btae119)
- [GMSC-mapper](https://github.com/BigDataBiology/GMSC-mapper) - [Python]

### databases
- [eggNOG 6.0](http://eggnog6.embl.de/) - [eggNOG 6.0: Enabling comparative genomics across 12,535 organisms.](https://doi.org/10.1093/nar/gkac1022)
- [KEGG](https://www.genome.jp/kegg/)
- [CARD](https://card.mcmaster.ca) - [CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database](https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/36263822/)
- [UniProt](https://www.uniprot.org) - [The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. _Nucleic Acids Res_. 51:D523–D531 (2023)](https://doi.org/10.1093/nar/gkac1052)

## microbial genome resources
- [GTDB](https://gtdb.ecogenomic.org/about) - [GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy](https://doi.org/10.1093/nar/gkab776)
- [GlobDB](https://globdb.org)
- [AMPSphere](https://ampsphere.big-data-biology.org/home) - [Discovery of antimicrobial peptides in the global microbiome with machine learning](https://doi.org/10.1016/j.cell.2024.05.013)
- [DRAMP](http://dramp.cpu-bioinfor.org) - [DRAMP 3.0: an enhanced comprehensive data repository of antimicrobial peptides](https://doi.org/10.1093/nar/gkab651)
- [proGenomes v3](https://progenomes.embl.de/index.cgi) - []
- [cFMD](https://github.com/SegataLab/cFMD) - [Carlino _et al_., Unexplored microbial diversity from 2,500 food metagenomes and links with the human microbiome, _Cell_, 2024](https://doi.org/10.1016/j.cell.2024.07.039)
- [GMSC](https://gmsc.big-data-biology.org) - [Duan, Y., Santos-Júnior, C.D., Schmidt, T.S. _et al_. A catalog of small proteins from the global microbiome. _Nat Commun_ 15, 7563 (2024)](https://doi.org/10.1038/s41467-024-51894-6)
- [GOMC](https://db.cngb.org/maya/datasets/MDB0000002) - [Chen, J., Jia, Y., Sun, Y. _et al_. Global marine microbial diversity and its potential in bioprospecting. _Nature_ 633, 371–379 (2024)](https://doi.org/10.1038/s41586-024-07891-2)
- [Ocean Microbiomics Database (OMD)](https://microbiomics.io/ocean/) - [Paoli, L., Ruscheweyh, HJ., Forneris, C.C. _et al_. Biosynthetic potential of the global ocean microbiome. _Nature_ **607**, 111–118 (2022)]( https://doi.org/10.1038/s41586-022-04862-3)
- [OceanDNA](https://doi.org/10.6084/m9.figshare.c.5564844.v1) - [Yoshizawa, Susumu; Nishimura, Yosuke (2022). The OceanDNA MAG catalog contains over 50,000 prokaryotic genomes originated from various marine environments. figshare. Collection](https://doi.org/10.6084/m9.figshare.c.5564844.v1)
- [AllTheBacteria](https://github.com/AllTheBacteria/AllTheBacteria) - [Martin Hunt, _et al_., _bioRxiv_, 2024. AllTheBacteria - all bacterial genomes assembled, available and searchable](https://doi.org/10.1101/2024.03.08.584059)
- [MGnify](https://www.ebi.ac.uk/metagenomics) - [Lorna Richardson, _et al_. MGnify: the microbiome sequence data analysis resource in 2023, _Nucleic Acids Research_ 51, D753-759 (2023)](https://doi.org/10.1093/nar/gkac1080]
- [Logan](https://registry.opendata.aws/pasteur-logan/)

## phylogenetics analysis
### phylogeny construction
- [IQ-TREE](http://www.iqtree.org) - [Bui Quang Minh, _et al_. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era, Molecular Biology and Evolution, Volume 37, Issue 5, May 2020, Pages 1530–1534](https://doi.org/10.1093/molbev/msaa015)
- [fasttree](http://www.microbesonline.org/fasttree/) - [Price, M.N., Dehal, P.S., and Arkin, A.P. (2010) FastTree 2 -- Approximately Maximum-Likelihood Trees for Large Alignments. _PLoS ONE_, 5(3):e9490](https://doi.org/10.1371/journal.pone.0009490)

### tree visualization
- [iTol](https://itol.embl.de) - [Web] - [Ivica Letunic, Peer Bork. Interactive Tree of Life (iTOL) v6: recent updates to the phylogenetic tree display and annotation tool](https://doi.org/10.1093/nar/gkae268)
- [ARB](http://www.arb-home.de) - [需要转发] - [Wolfgang Ludwig, _et al_. ARB: a software environment for sequence data. Nucleic Acids Research. 2004. 32(4):1363-1371](https://doi.org/10.1093/nar/gkh293)

## comparative genomics
- [ANI calculator](http://enve-omics.ce.gatech.edu/ani/) - [web]
- [AAI calculator](http://enve-omics.ce.gatech.edu/aai/) - [web]
- [fastANI](https://github.com/ParBLiSS/FastANI) - [C++] - [High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries](https://www.nature.com/articles/s41467-018-07641-9)
- [SVbyEye](https://github.com/daewoooo/SVbyEye) - [Davidd Porubsky, _et al_. SVbyEye: A visual tool to characterize structural variation among whole genome assemblies. _bioRxiv_ 2024](https://doi.org/10.1101/2024.09.11.612418)
- [CompareM](https://github.com/dparks1134/CompareM) - [Python] - Unsupported. Lasted version Dec 31 2020
- [EzAAI](https://github.com/endixk/ezaai) - [Java] - [Kim, D., Park, S. & Chun, J. Introducing EzAAI: a pipeline for high throughput calculations of prokaryotic average amino acid identity. _J Microbiol_. 59, 476–480 (2021)](https://doi.org/10.1007/s12275-021-1154-0)
- [SynTracker]( https://github.com/leylabmpi/SynTracker) - [Python] - [Enav, H., Paz, I. & Ley, R.E. Strain tracking in complex microbiomes using synteny analysis reveals per-species modes of evolution. _Nat Biotechnol_ (2024)](https://doi.org/10.1038/s41587-024-02276-2)

## statistical
- [LEfSe](https://github.com/SegataLab/lefse) - [Python] - [Nicola Segata, _et al_. Metagenomic Biomarker Discovery and Explanation. _Genome Biology_, 2011 Jun 24;12(6):R60](https://doi.org/10.1186%2Fgb-2011-12-6-r60)

## implemented workflows and platforms
- [Aviary](https://github.com/rhysnewell/aviary) - [Python] - a snakemake pipeline for assembly, binning, annotation, and strain diversity analysis.
- [MetaWRAP](https://github.com/bxlab/metaWRAP) - [Shell/Python] - [MetaWRAP—a flexible pipeline for genome-resolved metagenomic data analysis](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0541-1)
- [anvi'o](https://anvio.org) - [Eren, A.M., Kiefl, E., Shaiber, A. _et al_. Community-led, integrated, reproducible multi-omics with anvi’o. _Nat Microbiol_ 6, 3–6 (2021)](https://doi.org/10.1038/s41564-020-00834-3)
- [nano-rave](https://github.com/sanger-pathogens/nano-rave) - [Nextflow] - [Girgis, S.T., Adika, E., Nenyewodey, F.E. _et al_. Drug resistance and vaccine target surveillance of Plasmodium falciparum using nanopore sequencing in Ghana. _Nat Microbiol_ 8, 2365–2377 (2023)](https://doi.org/10.1038/s41564-023-01516-6)

## meta RNA
- [salmon](https://github.com/COMBINE-lab/salmon) - [Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods](https://doi.org/10.1038/nmeth.4197)

## Contact
This repository was created by [Jie Li](https://github.com/lijierr).

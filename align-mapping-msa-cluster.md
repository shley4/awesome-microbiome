## 比对 | alignment/mapping

### 短读长比对到长序列 | short-reads to sequence
- [[Bowtie2](https://github.com/BenLangmead/bowtie2)] - [C++, v2.5.4, 2024.05] - [2012.03, _Nat Methods_] - [[Fast gapped-read alignment with Bowtie 2.](https://doi.org/10.1038/nmeth.1923)]
- [[bwa](https://github.com/lh3/bwa)] - [C, v0.7.18, 2024.04] - [2009.05, _Bioinformatics_] - [[Fast and accurate short read alignment with Burrows-Wheeler transform.](https://doi.org/10.1093/bioinformatics/btp324)]
- [[bbmap](https://sourceforge.net/projects/bbmap/)] - [Java, v39.10, 2024.09] - [2014.03, _No. LBNL-7065E_] - [[BBMap: A Fast, Accurate, Splice-Aware Aligner](https://escholarship.org/uc/item/1h3515gn)]

### 长读长比对到长序列 | long-reads to sequence
- [BWA-MEM](https://github.com/lh3/bwa) - [C, v0.7.18, 2024.04] - [2013.5, _arXiv_] - [[Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM.](https://doi.org/10.48550/arXiv.1303.3997)]
- [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) - [C++, v2.2.1, 2021.3] - [2019.5, _IEEE Parallel and Distributed Processing Symposium (IPDPS)_] - [[Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems.](https://doi.org/10.1109/IPDPS.2019.00041)]
- [[Minimap2](https://github.com/lh3/minimap2)] - [C, v2.28, 2024.03] - [2018.05, _Bioinformatics_] - [[Minimap2: pairwise alignment for nucleotide sequences.](https://doi.org/10.1093/bioinformatics/bty191)]


### 序列比对到序列 | sequence to sequence
- [[blast+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)] - [C++, 2.16.0, 2024.06] - [2009, _BMC Bioinformatics_] - [[BLAST+: architecture and applications.](https://doi.org/10.1186/1471-2105-10-421)] - [[download BLAST databases.](https://ftp.ncbi.nlm.nih.gov/blast/db/)]
- [[DIAMOND](https://github.com/bbuchfink/diamond)] - [C++, v0.9.25及以下, 2019.07] - [2015.1, _Nat Methods_] - [[Fast and sensitive protein alignment using DIAMOND](https://doi.org/10.1038/nmeth.3176)]
- [[DIAMOND](https://github.com/bbuchfink/diamond)] - [C++, v2.1.9, 2024.01] - [2021.04, _Nat Methods_] - [[Sensitive protein alignments at tree-of-life scale using DIAMOND](https://doi.org/10.1038/s41592-021-01101-x)]
- [[LexicMap](https://github.com/shenwei356/LexicMap)] - [Go, v0.4.0, 2024.08] - [2024.08, _bioRxiv_] [[LexicMap: efficient sequence alignment against millions of prokaryotic genomes.](https://doi.org/10.1101/2024.08.30.610459)]
- [[foldmason](https://github.com/steineggerlab/foldmason)] - [C/C++, v1-763a428, 2024.08] - [2024.08, _bioRxiv_] - [[Multiple Protein Structure Alignment at Scale with FoldMason.](https://doi.org/10.1101/2024.08.01.606130)]

### 多序列比对 | multi sequence alignment
- [[MUSCLE](https://github.com/rcedgar/muscle)] - [C++, v5.2, 2024.08] - [2022.11, _Nat Commun_] - [[Muscle5: High-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny](https://doi.org/10.1038/s41467-022-34630-w)] - [Benchmark](https://github.com/rcedgar/muscle_benchmark)
- [[SSU-ALIGN](http://eddylab.org/software/ssu-align/)] - [Perl, v0.1.1, 2016.02] - [2009, _Ph.D thesis, Washington U_] - [[Structural RNA Homology Search and Alignment Using Covariance Models](http://eddylab.org/publications/Nawrocki09b/Nawrocki09b-phdthesis.pdf)]
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) - [C, v7.526, 2024.04] - [2013.01, _Mol Biol Evol_] - [[MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability](https://doi.org/10.1093/molbev/mst010)]
- [[Clustal 2](http://www.clustal.org/clustal2/)] - [NA, v2.1, 2010.10] - [2007.11, _Bioinformatics_] - [[Clustal W and Clustal X version 2.0](https://doi.org/10.1093/bioinformatics/btm404)]
- [[T-Coffee](https://github.com/cbcrg/tcoffee)] - [C/Perl, 13.46.1.b8b01e06, 2024.6] - [2000.09, _J Mol Biol_] - [[T-Coffee: a novel method for fast and accurate multiple sequence alignment](https://doi.org/10.1006/jmbi.2000.4042)]

### 序列聚类 | sequence cluster
- [MMseq](https://github.com/soedinglab/MMseqs) - [C++, Not Dev Anymore, 2016.9] - [2016.5, _Bioinformatics_] - [MMseqs software suite for fast and deep clustering and searching of large protein sequence sets](https://doi.org/10.1093/bioinformatics/btw006)
- [MMseq2](https://github.com/soedinglab/MMseqs2) - [C/C++, 15-6f452, 2023.10] - [2017.11, _Nat Biotechnol_] - [MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets](https://doi.org/10.1038/nbt.3988)
- [USEARCH](https://github.com/rcedgar/usearch12) - [C++/C, v12.0-beta1, 2024.6] - [2010.10, _Bioinformatics_] - [Search and clustering orders of magnitude faster than BLAST](https://doi.org/10.1093/bioinformatics/btq461)
- [USEARCH](http://www.drive5.com/usearch/) - [C++/C, usearch11.0.667, 看不出来] - [2010.10, _Bioinformatics_] - [Search and clustering orders of magnitude faster than BLAST](https://doi.org/10.1093/bioinformatics/btq461)
- [[CD-HIT](https://github.com/weizhongli/cdhit)] - [Perl/C++, v4.8.1, 2019.3] [2012.10, _Bioinformatics_] - [[CD-HIT: accelerated for clustering the next-generation sequencing data](https://doi.org/10.1093/bioinformatics/bts565)]

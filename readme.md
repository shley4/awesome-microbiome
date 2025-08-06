# awesome_microbiome
This repository serves as a continuously updated collection of algorithms, tools, databases, and tutorials for microbiome research.

## Table of content
- [Amplicon data analysis](#amplicon-data-analysis)
    - [Tools for amplicon](#tools-for-amplicon)
    - [rRNA Databases](#rrna-databases)
- [End-to-end workflows for metagenomes and isolates](#end-to-end-workflows-for-metagenomes-and-isolates)
    - [Workflows for metagenome](#workflows-for-metagenome)
    - [Workflows for isolate](#workflows-for-isolate)
- [Microbial genomes resources](#microbial-genomic-resources)
    - [Prokaryotic genome](#prokaryotic-genome)
    - [Phage genome](#phage-genome)
- [Quality control of sequencing data](#quality-control-of-sequencing-data)
    - [Reads simulation](#reads-simulation)
    - [Basecall](#basecall)
    - [QC](#qc)
    - [Correct reads](#correct-reads)
    - [Consensus sequence from long-reads](#consensus-sequence-from-long-reads)
- [Alignment and mapping](#alignment-and-mapping)
    - [Short-read to sequence](#short-read-to-sequence)
    - [Long-read to sequence](#long-read-to-sequence)
    - [Sequence to sequence](#sequence-to-sequence)
    - [Multi sequence alignment](#multi-sequence-alignment)
    - [Cluster sequence](#cluster-sequence)
- [Metagenome assmeble](#metagenome-assemble)
    - [Short read only assembly](#short-read-only-assembly)
    - [Long read only assmebly](#long-read-only-assembly)
    - [Hybrid assembly](#hybrid-assembly)
    - [PB HiFi read assembly](#pb-hifi-read-assembly)
    - [Polish](#polish)
    - [Genome size estimation](#genome-size-estimation)
    - [Other elements assembly](#other-elements-assmebly)
    - [Assembly improvement](#assembly-improvement)
- [Binning](#binning)
    - [Tools and workflows](#tools-and-workflows)
    - [Strain-level resolve](#strain-level-resolve)
    - [MAG improvement](#mag-improvement)
    - [MAG assessment](#mag-assessment)
    - [Other types of genome recovery](#other-genomes-recovery)
- [Gene and element prediction](#gene-and-element-prediction)
    - [ORF](#orf)
    - [Non-coding RNA](#non-coding-rna)
    - [Signal peptides](#signal-peptides)
    - [Plasmid databases](#plasmid-databases)
- [Taxonomy profile](#taxonomy-profile)
    - [Profile read](#profile-read)
    - [Profile contig](#profile-contig)
    - [MAG taxonomy](#mag-taxonomy)
- [Annotation](#annotation)
    - [Tools for annotation](#tools-for-annotation)
    - [Databases for annotation](#databases-for-annotation)
    - [Genome content analysis](#gennome-content-analysis)
    - [Protein structure analysis](#protein-structure-analysis)
    - [Elements](#elements)
    - [Bacteriophage](#bacteriophage)
- [Metabolic construction](#metabolic-construction)
    - [Tools for metabolic analysis](#tools-for-metabolic-analysis)
    - [Metabolic database](#metabolic-databases)
- [Comparative genomics](#comparative-genomics)
    - [AAI and ANI](#aai-and-ani)
    - [View comparative map](#view-comparative-map)
    - [HGT](#hgt)
    - [SV](#sv)
- [Visualization](#visualization)
    - [View MSA](#view-msa)
    - [View genome](#view-genome)
    - [View assemblies](#view-assemblies)
- [Phylogenetics](#phylogenetics)
    - [Build a tree](#build-a-tree)
    - [View tree](#view-tree)
- [Microbial diversity analysis](#microbial-diversity-analysis)
    - [Abundance](#abundance)
    - [Diversity](#diversity)
    - [Network](#network)
    - [Interactions](#interactions)
- [Metatranscriptomics](#metatranscriptomics)
    - [RNA quantification](#rna-quantification)
    - [RNA assembly](#rna-assembly)
    - [RNA SV](#rna-sv)
    - [Stats](#stats)
- [Surveillance](#surveillance)
- [Modifications](#modifications)
- [Pangenome related](#pangenome-related)


## Amplicon data analysis
### Tools for amplicon
- [`UPARSE`](http://drive5.com/uparse/) - (v12.0-beta1, 2024.6)
    - Edgar R. UPARSE: highly accurate OTU sequences from microbial amplicon reads. **_Nat Methods_** 10, 996–998 **(2013)**. https://doi.org/10.1038/nmeth.2604

- [`Mothur`](https://github.com/mothur/mothur/tree/v1.48.1) - (C++, v1.48.1, 2024.5)
    - Schloss PD, Westcott SL, Ryabin T, et al. Introducing mothur: Open-Source, Platform-Independent, Community-Supported Software for Describing and Comparing Microbial Communities. **_Appl Environ Microbiol_** 75 **(2009)**. https://doi.org/10.1128/AEM.01541-09

- [`dada2`](https://github.com/benjjneb/dada2) - (R, 1.26, 2022.11)
    - Callahan B, McMurdie P, Rosen M, et al. DADA2: High-resolution sample inference from Illumina amplicon data. **_Nat Methods_** 13, 581–583 **(2016)**. https://doi.org/10.1038/nmeth.3869

- [`QIIME2`](https://qiime2.org/) - (Java, 2024-05, 2024.5)
    - Bolyen E, Rideout JR, Dillon MR, et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. **_Nat Biotechnol_** 37, 852–857 **(2019)**. https://doi.org/10.1038/s41587-019-0209-9

- [`RDP Classifier`](https://sourceforge.net/projects/rdp-classifier/) - (Java, v2.14, 2023.8)
    - Wang Q, Cole JR. Updated RDP taxonomy and RDP Classifier for more accurate taxonomic classification. **_Microbiol Resour Announc_** 13, e01063-23 **(2024)**. https://doi.org/10.1128/mra.01063-23

- [`Tax4Fun1`](https://tax4fun.gobics.de) - [`Tax4Fun2`](https://github.com/fjossandon/Tax4Fun2) - (R, v1.1.6, 2019.11)
    - Aßhauer KP, Wemheuer B, Daniel R, Meinicke P. Tax4Fun: predicting functional profiles from metagenomic 16S rRNA data, **_Bioinformatics_** 31(17), 2882–2884 **(2015)**. https://doi.org/10.1093/bioinformatics/btv287
    - Wemheuer F, Taylor JA, Daniel R, et al. Tax4Fun2: prediction of habitat-specific functional profiles and functional redundancy based on 16S rRNA gene sequences. ***_Environmental Microbiome_*** 15(11) **(2020)**. https://doi.org/10.1186/s40793-020-00358-7

- [`PICRUSt`](https://github.com/picrust/picrust) - (Python, v1.1.4, 2019.6)
    - Langille M, Zaneveld J, Caporaso J, et al. Predictive functional profiling of microbial communities using 16S rRNA marker gene sequences. **_Nat Biotechnol_** 31, 814–821 **(2013)**. https://doi.org/10.1038/nbt.2676

- [`PICRUSt2`](https://github.com/picrust/picrust2) - (Python, v2.6.2, 2025.4)
    - Douglas GM, Maffei VJ, Zaneveld JR, et al. PICRUSt2 for prediction of metagenome functions. **_Nat Biotechnol_** 38, 685–688 **(2020)**. https://doi.org/10.1038/s41587-020-0548-6

- [`graftM`](https://github.com/geronimp/graftM) - (Python, v0.14.0, 2022.5)
    - Boyd JA, Woodcroft BJ, Tyson GW. GraftM: a tool for scalable, phylogenetically informed classification of genes within metagenomes. **_Nucleic Acids Res_** 46(10), e59 **(2018)**. https://doi.org/10.1093/nar/gky174

- [`RiboSnake`](https://github.com/IKIM-Essen/RiboSnake) - (snakemake, v0.10.0, 2024.8)
    - Dörr AK, Welling J, Dörr A, et al. RiboSnake – a user-friendly, robust, reproducible, multipurpose and documentation-extensive pipeline for 16S rRNA gene microbiome analysis. **_Gigabyte_** **(2024)**. https://doi.org/10.46471/gigabyte.132

- [`ssUMI`](https://github.com/ZielsLab/ssUMI) - (Shell, NoReleaseTag)
    - Lin X, Waring K, Ghezzi H, et al. High accuracy meets high throughput for near full-length 16S ribosomal RNA amplicon sequencing on the Nanopore platform. **_PNAS Nexus_** 3(10), pgae411 **(2024)**. https://doi.org/10.1093/pnasnexus/pgae411


### rRNA Databases
- [`EUKARYOME`](https://eukaryome.org) - (18S/ITS/28S, v1.9.2, 2024.8)
    - Tedersoo L, Moghaddam MSH, Mikryukov V, et al. EUKARYOME: the rRNA gene reference database for identification of all eukaryotes, **_Database_** 2024, baae043 **(2024)**. https://doi.org/10.1093/database/baae043

- [`Silva`](https://www.arb-silva.de) - (rRNA, v138.2, 2024.7) - A continuously updated rRNA database with high update frequency, currently the largest and most comprehensive of its kind.
    - Quast C, Pruesse E, Yilmaz P, et al. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. **_Nucleic Acids Res_** 41(D1), D590–D596 **(2013)**. https://doi.org/10.1093/nar/gks1219

- [`Greengenes2`](http://ftp.microbio.me/greengenes_release/) - (16S, v2022.10, 2022.10) - The first version was last updated in 2012. Greengenes2 was released in 2022 and is updated approximately every two years. 
    - DeSantis TZ, Hugenholtz P, Larsen N, et al. Greengenes, a chimera-checked 16S rRNA gene database and workbench compatible with ARB. **_Appl Environ Microbiol_** 72(7), 5069-72 **(2006)**. https://doi.org/10.1128/AEM.03006-05
    - McDonald D, Jiang Y, Balaban M, et al. Greengenes2 unifies microbial data in a single reference tree. **_Nat Biotechnol_** 42, 715–718 **(2024)**. https://doi.org/10.1038/s41587-023-01845-1

- [`MiDAS`](https://midasfieldguide.org/guide) - (16S, v5, 2024.7)
    - Dueholm MKD, Andersen KS, Korntved AKC, et al. MiDAS 5: Global diversity of bacteria and archaea in anaerobic digesters. **_Nat Commun_** 15, 5361 **(2024)**. https://doi.org/10.1038/s41467-024-49641-y

- [`KSGP`](https://ksgp.earlham.ac.uk/index.php?site=download) - (SSU, v3.1)
    - Grant A, Aleidan A, Davies CS, et al. KSGP 3.1: improved taxonomic annotation of Archaea communities using LotuS2, the Genome Taxonomy Database and RNAseq data. **_ISME Communications_** ycaf094 **(2025)**. https://doi.org/10.1093/ismeco/ycaf094

## End-to-end workflows for metagenome and isolate
### Workflows for metagenome
- [`Aviary`](https://github.com/rhysnewell/aviary) - (Snakemake, v0.9.2, 2024.09)
    - A pipeline for assembly, binning, annotation, and strain diversity analysis.

- [`MAG`](https://github.com/nf-core/mag) - (Nextflow, v3.0.3, 2024.8)
    - Krakau S, Straub D, Gourlé H, et al. nf-core/mag: a best-practice pipeline for metagenome hybrid assembly and binning. **_NAR Genomics and Bioinformatics_** 4(1), lqac007 **(2022)**. https://doi.org/10.1093/nargab/lqac007

- [`MetaWRAP`](https://github.com/bxlab/metaWRAP) - (Shell/Python, v1.3, 2020.08) - A modular pipeline for metagenomic analysis, covering steps such as quality control, assembly, binning, genome refinement, classification, and annotation. Users can independently run any individual module.
    - Uritskiy GV, DiRuggiero J & Taylor J. MetaWRAP—a flexible pipeline for genome-resolved metagenomic data analysis. **_Microbiome_** 6, 158 **(2018)**. https://doi.org/10.1186/s40168-018-0541-1

- [`VEBA`](https://github.com/jolespin/veba) - (Python, v2.3.0, 2024.9)
    - Espinoza JL, Dupont CL. VEBA: a modular end-to-end suite for in silico recovery, clustering, and analysis of prokaryotic, microeukaryotic, and viral genomes from metagenomes. **_BMC Bioinformatics_** 23, 419 **(2022)**. https://doi.org/10.1186/s12859-022-04973-8

- [`ATLAS`](https://github.com/metagenome-atlas/atlas) - (Python, v2.19.0, 2024.7)
    - Kieser S, Brown J, Zdobnov EM, et al. ATLAS: a Snakemake workflow for assembly, annotation, and genomic binning of metagenome sequence data. **_BMC Bioinformatics_** 21, 257 **(2020)**. https://doi.org/10.1186/s12859-020-03585-4

- [`anvi'o`](https://github.com/merenlab/anvio/tree/v8) - (Win/MacOS/Linux, v8, 2023.09)
    - Eren AM, Kiefl E, Shaiber A, et al. Community-led, integrated, reproducible multi-omics with anvi’o. **_Nat Microbiol_** 6, 3–6 **(2021)**. https://doi.org/10.1038/s41564-020-00834-3

- [`nano-rave`](https://github.com/sanger-pathogens/nano-rave) - (Nextflow, v1.0.0, 2023.01)
    - Girgis ST, Adika E, Nenyewodey FE, et al. Drug resistance and vaccine target surveillance of Plasmodium falciparum using nanopore sequencing in Ghana. **_Nat Microbiol_** 8, 2365–2377 **(2023)**. https://doi.org/10.1038/s41564-023-01516-6

- [`metagWGS`](https://forgemia.inra.fr/genotoul-bioinfo/metagwgs) - (Nextflow, v2.4.2, 2023.6)
    - Mainguy J, Vienne M, Fourquet J, et al. metagWGS, a comprehensive workflow to analyze metagenomic data using Illumina or PacBio HiFi reads. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.09.13.612854

- [`SqueezeMeta`](https://github.com/jtamames/SqueezeMeta) - (C/C++/Python/Perl, v1.6.5, 2024.8)
    - Tamames J, Puente-Sánchez F. SqueezeMeta, A Highly Portable, Fully Automatic Metagenomic Analysis Pipeline. **_Front Microbiol_** 9 **(2019)**. https://doi.org/10.3389/fmicb.2018.03349

- [`slamM`](https://github.com/Ecogenomics/slamM/tree/master) - (Snakemake)

- [`CAT_pack`](https://github.com/MGXlab/CAT_pack) - (Python, v6.0.1, 2024.3)
    - von Meijenfeldt FAB, Arkhipova K, Cambuy DD, et al. Robust taxonomic classification of uncharted microbial sequences and bins with CAT and BAT. **_Genome Biol_** 20, 217 **(2019)**. https://doi.org/10.1186/s13059-019-1817-x

- [`Metagenomics-Toolkit`](https://github.com/metagenomics/metagenomics-tk) - (Nextflow, v0.4.5, 2024.10)
    - Belmann P, Osterholz B, Kleinbölting N, et al. Metagenomics-Toolkit: The Flexible and Efficient Cloud-Based Metagenomics Workflow featuring Machine Learning-Enabled Resource Allocation. bioRxiv **(2024)**. https://doi.org/10.1101/2024.10.22.619569

- [`BugBuster`](https://github.com/gene2dis/BugBuster) - Nextflow
    - Fuentes-Santander F, Curiqueo C, Araos R, Ugalde JA. BugBuster: A novel automatic and reproducible workflow for metagenomic data analysis. **bioRxiv** **(2025)**. https://doi.org/10.1101/2025.02.24.639915

- [`MARTi`](https://marti.cyverseuk.org/)
    - Peel N, Martin S, Heavens D, et al. MARTi: a real-time analysis and visualisation tool for nanopore metagenomics. **bioRxiv** **(2025)**. https://doi.org/10.1101/2025.02.14.638261

### Workflows for isolate
- [`Bactopia`](https://github.com/bactopia/bactopia) - (Nextflow/Perl, v3.1.0, 2024.9)
    - Petit RA and Read TD. Bactopia: a Flexible Pipeline for Complete Analysis of Bacterial Genomes. **_mSystems_** 5 **(2020)**. https://doi.org/10.1128/msystems.00190-20

- [`ASA³P`](https://github.com/oschwengers/asap) - (Groovy/JS, v1.3.0, 2020.5)
    - Schwengers O, Hoek A, Fritzenwanker M, et al. ASA3P: An automatic and scalable pipeline for the assembly, annotation and higher-level analysis of closely related bacterial isolates. **_PLoS Comput Biol_** 16(3), e1007134 **(2020)**. https://doi.org/10.1371/journal.pcbi.1007134

- [`microPIPE`](https://github.com/BeatsonLab-MicrobialGenomics/micropipe) - (html tutorial Nextflow, step by step)
    - Murigneux V, Roberts LW, Forde BM, et al. MicroPIPE: validating an end-to-end workflow for high-quality complete bacterial genome construction. **_BMC Genomics_** 22, 474 **(2021)**. https://doi.org/10.1186/s12864-021-07767-z

- [`AQUAMIS`](https://gitlab.com/bfr_bioinformatics/AQUAMIS) - (Python/Shell, v1.4.2, 2024.6)
    - Deneke C, Brendebach H, Uelze L, et al. Species-Specific Quality Control, Assembly and Contamination Detection in Microbial Isolate Sequences with AQUAMIS. **_Genes_** 12(5), 644 **(2021)**. https://doi.org/10.3390/genes12050644

- [`Nullarbor`](https://github.com/tseemann/nullarbor) - (Perl, v1.41, 2018.6) - Pipeline to generate complete public health microbiology reports from sequenced isolates.

- [`ProkEvo`](https://github.com/npavlovikj/ProkEvo) - jupyter tutorial
    - Pavlovikj N, Gomes-Neto JC, Deogun JS, Benson AK. ProkEvo: an automated, reproducible, and scalable framework for high-throughput bacterial population genomics analyses. **_PeerJ_** 9, e11376 **(2021)**. https://doi.org/10.7717/peerj.11376

- [`Public Health Bioinformatics`](https://github.com/theiagen/public_health_bioinformatics) - (WDL, v3.1.0, 2025.7) - TheiaCoV(viral), TheiaProk(bacterial)
    - Libuit KG, Doughty EL, Otieno JR, et al. Accelerating Bioinformatics Implementation in Public Health. **_Microbial Genomics_** 9(7), 2057-5858 **(2023)**. https://doi.org/10.1099/mgen.0.001051

- [`Public Health Bioinformatics`](https://github.com/theiagen/public_health_bioinformatics) - (WDL, v3.1.0, 2025.7) - TheiaEuk
    - Ambrosio FJ, Scribner MR, Wright SM, et al. TheiaEuk: A Species-Agnostic Bioinformatics Workflow for Fungal Genomic Characterization. **_Frontiers in Public Health_** 11 **(2023)**. https://doi.org/10.3389/fpubh.2023.1198213

- [`rMAP`](https://github.com/GunzIvan28/rMAP) - (Shell/Python, tutorial)
    - Sserwadda I, Mboowa G. rMAP: the Rapid Microbial Analysis Pipeline for ESKAPE bacterial group whole-genome sequence data. **_Microbial Genomics_** 7 **(2021)**. https://doi.org/10.1099/mgen.0.000583

- [`TORMES`](https://github.com/nmquijada/tormes) - (Shell, v1.3.0, 2021.8)
    - Quijada NM, Rodríguez-Lázaro D, Eiros JM, Hernández M. TORMES: an automated pipeline for whole bacterial genome analysis. **_Bioinformatics_** 35(21), 4207–4212 **(2019)**. https://doi.org/10.1093/bioinformatics/btz220

- [`Hybracter`](https://github.com/gbouras13/hybracter) - (Python, v0.10.0, 2024.10)
    - Bouras G, Houtak G, Wick R, et al. Hybracter: Enabling Scalable, Automated, Complete and Accurate Bacterial Genome Assemblies. **_Microbial Genomics_** 10, 5 **(2024)**. https://doi.org/10.1099/mgen.0.001244

## Microbial genomic resources
### Prokaryotic genome
- [`GTDB`](https://gtdb.ecogenomic.org/about) - (Genomes, R220, 2024.04)
    - Parks DH, Chuvochina M, Rinke C, et al. GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. **_Nucleic Acids Res_** 50(D1), D785–D794 **(2022)**. https://doi.org/10.1093/nar/gkab776

- [`GlobDB`](https://globdb.org) - (Genomes, 220) - A collection from Daan Septh
    - Speth DR, Pullen N, Aroney STN, et al. GlobDB: A comprehensive species-dereplicated microbial genome resource. **arXiv** **(2025)**. https://doi.org/10.48550/arXiv.2506.11896

- [`proGenomes`](https://progenomes.embl.de/index.cgi) - (Genomes, v3)

- [`cFMD`](https://github.com/SegataLab/cFMD) - (Food Microb Genomes, v1.2.0, 2024.08)
    - Carlino N, Blanco-Míguez A, Punčochář M, et al. Unexplored microbial diversity from 2,500 food metagenomes and links with the human microbiome. **_Cell_** 187(20), 5775-5795 **(2024)**. https://doi.org/10.1016/j.cell.2024.07.039

- [`GOMC`](https://db.cngb.org/maya/datasets/MDB0000002) - (Genomes, 2024.09)
    - Chen J, Jia Y, Sun Y, et al. Global marine microbial diversity and its potential in bioprospecting. **_Nature_** 633, 371–379 **(2024)**. https://doi.org/10.1038/s41586-024-07891-2

- [`Ocean Microbiomics Database (OMD)`](https://microbiomics.io/ocean/) - (Genomes, v1.1)
    - Paoli L, Ruscheweyh HJ, Forneris CC, et al. Biosynthetic potential of the global ocean microbiome. **_Nature_** 607, 111–118 **(2022)**. https://doi.org/10.1038/s41586-022-04862-3

- [`OceanDNA`](https://doi.org/10.6084/m9.figshare.c.5564844.v1) - (Genomes, v2022-05, 2022.05)
    - Yoshizawa S, Nishimura Y. The OceanDNA MAG catalog contains over 50,000 prokaryotic genomes originated from various marine environments. **_figshare_** Collection **(2022)**. https://doi.org/10.6084/m9.figshare.c.5564844.v1

- [`AllTheBacteria`](https://github.com/AllTheBacteria/AllTheBacteria) - (Genomes, v2024-08, 2024.08)
    - Hunt M, Lima L, Anderson D, et al. AllTheBacteria - all bacterial genomes assembled, available and searchable. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.03.08.584059

- [`MGnify`](https://www.ebi.ac.uk/metagenomics) - (Genomes, 2023.09) - ftp: https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/
    - Richardson L, Allen B, Baldi G, et al. MGnify: the microbiome sequence data analysis resource in 2023. **_Nucleic Acids Res_** 51(D1), D753–D759 **(2023)**. https://doi.org/10.1093/nar/gkac1080

- [`MGnify genome`](https://www.ebi.ac.uk/metagenomics/browse/genomes/)
    - Gurbich TA, Almeida A, Beracochea M, et al. MGnify Genomes: A Resource for Biome-specific Microbial Genome Catalogues. **_J Mol Biol_** 435, 14 **(2023)**. https://doi.org/10.1016/j.jmb.2023.168016

- [`Logan`](https://registry.opendata.aws/pasteur-logan/) - (Unitigs/Contigs, AWS)

- [`BakRep`](https://bakrep.computational.bio)
    - Fenske L, Jelonek L, Goesmann A, Schwengers O. BakRep – a searchable large-scale web repository for bacterial genomes, characterizations and metadata. **_Microbial Genomics_** 10, 10 **(2024)**. https://doi.org/10.1099/mgen.0.001305

- [`BacDive`](https://bacdive.dsmz.de/) - Web: https://bacdive.dsmz.de
    - Schober I, Koblitz J, Carbasse JS, et al. BacDive in 2025: the core database for prokaryotic strain data. **_Nucleic Acids Res_** 53(D1), D748–D756 **(2025)**. https://doi.org/10.1093/nar/gkae959

- [`GROWdb, NMDC data portal`](https://data.microbiomedata.org/)/[with KBase](https://narrative.kbase.us/collections/GROW)/[GROWdb Explorer](https://geocentroid.shinyapps.io/GROWdatabase/)
    - Borton MA, McGivern BB, Willi KR, et al. A functional microbiome catalogue crowdsourced from North American rivers. **_Nature_** 637, 103–112 **(2025)**. https://doi.org/10.1038/s41586-024-08240-z

### Phage genome
- [`inphared`](https://github.com/RyanCook94/inphared)
    - Cook R, Brown N, Redgwell T, et al. INfrastructure for a PHAge REference Database: Identification of Large-Scale Biases in the Current Collection of Cultured Phage Genomes. **_Phage_** 2, 4 **(2021)**. https://doi.org/10.1089/phage.2021.0007

- [`PhageDive`](https://phagedive.dsmz.de/) - Web Search
    - Rolland C, Wittmann J, Reimer LC, et al. PhageDive: the comprehensive strain database of prokaryotic viral diversity, **_Nucleic Acids Res_** 53(D1), D819–D825 **(2025)**. https://doi.org/10.1093/nar/gkae878

## Quality control of sequencing data
### Reads simulation
- [`wgsim`](https://github.com/lh3/wgsim) - (C) Reads simulator

- [`badread`](https://github.com/rrwick/Badread) - (Python, v0.4.1, 2024.2) - a long read simulator that can imitate many types of read problems

### Basecall
- [`Dorado`](https://github.com/nanoporetech/dorado) - (C++, v0.9.0, 2024.12) - Oxford Nanopore's Basecaller

### QC
- [`FastQC`](https://github.com/s-andrews/FastQC) - (Java, 0.12.1, 2023.3) - A quality control tool for high throughput sequence data

- [`fastp`](https://github.com/OpenGene/fastp) - (C++/C, v0.23.4, 2023.5)
    - Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. **_Bioinformatics_** 34(17), 884–890 **(2018)**. https://doi.org/10.1093/bioinformatics/bty560
    - Chen S. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. **_iMeta_** 2(2), e107 **(2023)**. https://doi.org/10.1002/imt2.107

- [`Trimmomatic`](https://github.com/usadellab/Trimmomatic) - (Java, v0.39, 2021.01)
    - Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. **_Bioinformatics_** 30(15), 2114–2120 **(2014)**. https://doi.org/10.1093/bioinformatics/btu170

- [`Cutadapt`](https://github.com/marcelm/cutadapt) - (Python, v4.9, 2024.6)
    - Martin M. Cutadapt Removes Adapter Sequences From High-Throughput Sequencing Reads. **_EMBnet J_** 17, 1 **(2010)**. https://doi.org/10.14806/ej.17.1.200

- [`FastUniq`](https://sourceforge.net/projects/fastuniq/files/) - (C++, 1.1, 2012.05)
    - Xu H, Luo X, Qian J, et al. FastUniq: A Fast De Novo Duplicates Removal Tool for Paired Short Reads. **_PLoS ONE_** 7(12), e52249 **(2012)**. https://doi.org/10.1371/journal.pone.0052249

- [`PoreChop`](https://github.com/rrwick/Porechop) - (C++, v0.2.4, 2018.10) - adapter trimmer for Oxford Nanopore reads. Unsupported since Oct 2018

- [`Sickle`](https://github.com/najoshi/sickle) - (C, v1.33, 2014.7) - Sickle: A sliding-window, adaptive, quality-based trimming tool for FastQ files

- [`PEAT`](https://github.com/jhhung/PEAT) - (C++, v1.2.4, 2016.2) - No longer maintained, recommond EARRINGS 
    - Li YL, Weng JC, Hsiao CC. et al. PEAT: an intelligent and efficient paired-end sequencing adapter trimming algorithm.**_BMC Bioinformatics_** 16(Suppl 1), S2 **(2015)**. https://doi.org/10.1186/1471-2105-16-S1-S2

- [`EARRINGS`](https://github.com/jhhung/EARRINGS) - (C++, v1.1.0, 2024.6)
    - Wang TH, Huang CC, Hung JH. EARRINGS: an efficient and accurate adapter trimmer entails no a priori adapter sequences. **_Bioinformatics_** 37(13), 1846–1852 **(2021)**. https://doi.org/10.1093/bioinformatics/btab025

- [`NanoPack`](https://github.com/wdecoster/nanopack) - (Python) - No release
    - De Coster W, D’Hert S, Schultz DT, et al. NanoPack: visualizing and processing long-read sequencing data. **_Bioinformatics_** 34(15), 2666–2669 **(2018)**. https://doi.org/10.1093/bioinformatics/bty149

- [`NanoPlot`](https://github.com/wdecoster/NanoPlot) - (Python, v1.42.0, 2023.10)
    -  De Coster W, Rademakers R. NanoPack2: population-scale evaluation of long-read sequencing data. **_Bioinformatics_** 39(5), btad311 **(2023)**. https://doi.org/10.1093/bioinformatics/btad311

- [`chopper`](https://github.com/wdecoster/chopper) - (Rust, v0.9.0, 2024.8) - cite NanoPack2
    -  De Coster W, Rademakers R. NanoPack2: population-scale evaluation of long-read sequencing data. **_Bioinformatics_** 39(5), btad311 **(2023)**. https://doi.org/10.1093/bioinformatics/btad311

- [`nanoQC`](https://github.com/wdecoster/nanoQC) - (Python) - no release
    - De Coster W, D’Hert S, Schultz DT, et al. NanoPack: visualizing and processing long-read sequencing data. **_Bioinformatics_** 34(15), 2666–2669 **(2018)**. https://doi.org/10.1093/bioinformatics/bty149

- [`fastplong`](https://github.com/OpenGene/fastplong) - (C++, v0.2.2, 2024.12) - Ultra-fast preprocessing and quality control for long-read sequencing data

- [`sequali`](https://github.com/rhpvorderman/sequali) - (C/python, v0.12.0, 2024.10)
    - Vorderman RHP. Sequali: efficient and comprehensive quality control of short- and long-read sequencing data. **_Bioinformatics Advances_** 5(1), vbaf010 **(2025)**. https://doi.org/10.1093/bioadv/vbaf010

### Correct reads
### Consensus sequence from long-reads
- [`medaka`](https://github.com/nanoporetech/medaka) - (Python, v2.0.0, 2024.9) - a tool to create consensus sequences and variant calls from nanopore sequencing data

## Alignment and mapping
### Short-read to sequence
- [`Bowtie2`](https://github.com/BenLangmead/bowtie2) - (C++, v2.5.4, 2024.05)
    - Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. **_Nat Methods_** 9, 357–359 **(2012)**. https://doi.org/10.1038/nmeth.1923

- [`bwa`](https://github.com/lh3/bwa) - (C, v0.7.18, 2024.04)
    - Li H, Durbin R. Fast and accurate short read alignment with Burrows–Wheeler transform. **_Bioinformatics_** 25(14), 1754–1760 **(2009)**. https://doi.org/10.1093/bioinformatics/btp324

- [`bbmap`](https://sourceforge.net/projects/bbmap/) - (Java, v39.10, 2024.09)
    - Bushnell B. BBMap: A Fast, Accurate, Splice-Aware Aligner. **_LBNL Report_** **(2014)**. https://escholarship.org/uc/item/1h3515gn

- [`strobealign`](https://github.com/ksahlin/strobealign) - (C++, v0.15.0, 2024.12)
    - Sahlin K. Strobealign: flexible seed size enables ultra-fast and accurate read alignment. **_Genome Biol_** 23, 260 **(2022)**. https://doi.org/10.1186/s13059-022-02831-7

### Long-read to sequence
- [`BWA-MEM`](https://github.com/lh3/bwa) - (C, v0.7.18, 2024.04)
    - Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. **arXiv** **(2013)**. https://doi.org/10.48550/arXiv.1303.3997

- [`BWA-MEM2`](https://github.com/bwa-mem2/bwa-mem2) - (C++, v2.2.1, 2021.3)
    - Vasimuddin M, Misra S, Li H and Aluru S. Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. **_IEEE International Parallel and Distributed Processing Symposium (IPDPS)_** 314-324 **(2019)**. https://doi.org/10.1109/IPDPS.2019.00041

- [`Minimap2`](https://github.com/lh3/minimap2) - (C, v2.28, 2024.03)
    - Li H. Minimap2: pairwise alignment for nucleotide sequences. **_Bioinformatics_** 34(18), 3094–3100 **(2018)**. https://doi.org/10.1093/bioinformatics/bty191

### Sequence to sequence
- [`blast+`](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) - (C++, 2.16.0, 2024.06) - download databases: https://ftp.ncbi.nlm.nih.gov/blast/db/
    - Camacho C, Coulouris G, Avagyan V, et al. BLAST+: architecture and applications. **_BMC Bioinformatics_** 10, 421 **(2009)**. https://doi.org/10.1186/1471-2105-10-421

- [`DIAMOND`](https://github.com/bbuchfink/diamond) - (C++, v2.1.9, 2024.01)
    - Buchfink B, Xie C & Huson D. Fast and sensitive protein alignment using DIAMOND. **_Nat Methods_** 12, 59–60 **(2015)**. https://doi.org/10.1038/nmeth.3176
    - Buchfink B, Reuter K & Drost HG. Sensitive protein alignments at tree-of-life scale using DIAMOND. **_Nat Methods_** 18, 366–368 **(2021)**. https://doi.org/10.1038/s41592-021-01101-x

- [`LexicMap`](https://github.com/shenwei356/LexicMap) - (Go, v0.4.0, 2024.08)
    - Shen W, Iqbal Z. LexicMap: efficient sequence alignment against millions of prokaryotic genomes. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.08.30.610459

- [`foldmason`](https://github.com/steineggerlab/foldmason) - (C/C++, v1-763a428, 2024.08)
    - Gilchrist CLM, Mirdita M, Steinegger M. Multiple Protein Structure Alignment at Scale with FoldMason. **bioRixv** **(2024)**. https://doi.org/10.1101/2024.08.01.606130

- [`SkyBLAST`](https://sky-blast.com/) - a replica of the NIH's interface - providing an alternative to the US gov service that's less congested, faster & more reliable

### Multi sequence alignment
- [`MUSCLE`](https://github.com/rcedgar/muscle) - (C++, v5.2, 2024.08) - Benchmark:https://github.com/rcedgar/muscle_benchmark
    - Edgar RC. Muscle5: High-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny. **_Nat Commun_** 13, 6968 **(2022)**. https://doi.org/10.1038/s41467-022-34630-w

- [`SSU-ALIGN`](http://eddylab.org/software/ssu-align/) - (Perl, v0.1.1, 2016.02)
    - Nawrocki EP. Structural RNA Homology Search and Alignment Using Covariance Models. _Ph.D thesis, Washington U_ **(2009)**. http://eddylab.org/publications/Nawrocki09b/Nawrocki09b-phdthesis.pdf

- [`MAFFT`](https://mafft.cbrc.jp/alignment/software/) - (C, v7.526, 2024.04)
    - Katoh K, Standley DM. MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability. **_Mol Biol Evol_** 30(4), 772–780 **(2013)**. https://doi.org/10.1093/molbev/mst010

- [`Clustal 2`](http://www.clustal.org/clustal2/) - (v2.1, 2010.10)
    - Larkin MA, Blackshields G, Brown NP, et al. Clustal W and Clustal X version 2.0. **_Bioinformatics_** 23(21), 2947–2948 **(2007)**. https://doi.org/10.1093/bioinformatics/btm404

- [`T-Coffee`](https://github.com/cbcrg/tcoffee) - (C/Perl, 13.46.1.b8b01e06, 2024.6)
    - Notredame C, Higgins DG, Heringa J. T-Coffee: a novel method for fast and accurate multiple sequence alignment. **_J Mol Biol_** 302(1), 205-217 **(2000)**. https://doi.org/10.1006/jmbi.2000.4042

### Cluster sequence
- [`MMseq`](https://github.com/soedinglab/MMseqs) - (C++, 2016.9)
    - Hauser M, Steinegger M, Söding J. MMseqs software suite for fast and deep clustering and searching of large protein sequence sets. **_Bioinformatics_** 32(9), 1323–1330 **(2016)**. https://doi.org/10.1093/bioinformatics/btw006

- [`MMseq2`](https://github.com/soedinglab/MMseqs2) - (C/C++, 15-6f452, 2023.10)
    - Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. **_Nat Biotechnol_** 35, 1026–1028 **(2017)**. https://doi.org/10.1038/nbt.3988
    
- [`MMseq2 CUDA GPU-based`](https://github.com/soedinglab/mmseqs2) - NVIDIA blog: https://developer.nvidia.com/blog/boost-alphafold2-protein-structure-prediction-with-gpu-accelerated-mmseqs2
    - Kallenborn F, Chacon A, Hundt C, et al. GPU-accelerated homology search with MMseqs2. **bioRxiv** **(2025)**. https://doi.org/10.1101/2024.11.13.623350

- [`USEARCH`](https://github.com/rcedgar/usearch12) - (C++/C, v12.0-beta1, 2024.6)
    - Edgar RC. Search and clustering orders of magnitude faster than BLAST. **_Bioinformatics_** 26(19), 2460–2461 **(2010)**. https://doi.org/10.1093/bioinformatics/btq461

- [`CD-HIT`](https://github.com/weizhongli/cdhit) - (Perl/C++, v4.8.1, 2019.3)
    - Fu L, Niu B, Zhu Z. CD-HIT: accelerated for clustering the next-generation sequencing data. **_Bioinformatics_** 28(23), 3150–3152 **(2012)**. https://doi.org/10.1093/bioinformatics/bts565


## Metagenome assembly
### Short read only assembly
- [`SOAPdenovo2`](https://github.com/aquaskyline/SOAPdenovo2) - (C, r242, 2020.10)
    - Luo R, Liu B, Xie Y, et al. SOAPdenovo2: an empirically improved memory-efficient short-read de novo assembler. **_GigaScience_** 1(1), 2047–217X–1–18 **(2012)**. https://doi.org/10.1186/2047-217X-1-18

- [`idba/IDBA-UD`](https://github.com/loneknightpy/idba) - (C++, v1.1.3, 2016.7)
    - Peng Y, Leung HCM, Yiu SM, Chin FYL. IDBA-UD: a de novo assembler for single-cell and metagenomic sequencing data with highly uneven depth. **_Bioinformatics_** 28(11), 1420–1428 **(2012)**. https://doi.org/10.1093/bioinformatics/bts174

- [`SPAdes/metaSPAdes`](https://github.com/ablab/spades) - (C++, v4.0.0, 2024.6)
    - Sergey Nurk, Dmitry Meleshko, Anton Korobeynikov and Pavel A. Pevzner. metaSPAdes: a new versatile metagenomic assembler. **_Genome Res_** 27, 824-834 **(2017)**. https://doi.org/10.1101/gr.213959.116

- [`megahit`](https://github.com/voutcn/megahit) - (C++, v1.2.9, 2019.10) - A metagenome assembler with high speed and low memory consumption.
    - Li D, Liu CM, Luo R, et al. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph **_Bioinformatics_** 31(10), 1674–1676 **(2015)**. https://doi.org/10.1093/bioinformatics/btv033
    - Li D, Luo R, Liu CM, et al. MEGAHIT v1.0: A Fast and Scalable Metagenome Assembler driven by Advanced Methodologies and Community Practices. **_Methods_** 102, 3-11 **(2016)**. https://doi.org/10.1016/j.ymeth.2016.02.020

- [`PenguiN`](https://github.com/soedinglab/plass) - (C, v5-cf8933, 2024.3)
    - Jochheim A, Jochheim FA, Kolodyazhnaya A, et al. Strain-resolved de-novo metagenomic assembly of viral genomes and microbial 16S rRNAs. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.03.29.587318

- [`BBMerge`](https://sourceforge.net/projects/bbmap/) - (Bash, v39.10, 2024.9)
    - Bushnell B, Rood J, Singer E. BBMerge – Accurate paired shotgun read merging via overlap. **_PLoS ONE_** 12(10), e0185056 **(2017)**. https://doi.org/10.1371/journal.pone.0185056

- [`FLASH`](https://github.com/ebiggers/flash) - (C, v1.2.11, 2018.3) 
    - Magoč T, Salzberg SL. FLASH: fast length adjustment of short reads to improve genome assemblies. **_Bioinformatics_** 27(21), 2957–2963 **(2011)**. https://doi.org/10.1093/bioinformatics/btr507

- [`ABySS`](https://github.com/bcgsc/abyss) - (C++/C, v2.3.9, 2024.9)
    - Jackman SD, Vandervalk BP, Mohamadi H, et al. ABySS 2.0: resource-efficient assembly of large genomes using a Bloom filter. **_Genome Res_** 27, 768-777 **(2017)**. http://www.genome.org/cgi/doi/10.1101/gr.214346.116

### Long read only assembly
- [`MaSuRCA`](https://github.com/alekseyzimin/masurca) - (M4/Perl, v4.1.1, 2024.2) 
    - Zimin AV, Marçais G, Puiu D, et al. The MaSuRCA genome assembler. **_Bioinformatics_** 29(21), 2669–2677 **(2013)**. https://doi.org/10.1093/bioinformatics/btt476
    - Zimin AV, Salzberg SL. The SAMBA tool uses long reads to improve the contiguity of genome assemblies. **_PLoS Comput Biol_** 18(2), e1009860 **(2022)**. https://doi.org/10.1371/journal.pcbi.1009860

- [`metaFlye`](https://github.com/fenderglass/Flye) - (C/C++, v2.9.5, 2024.8)
    - Kolmogorov M, Bickhart DM, Behsaz B, et al. metaFlye: scalable long-read metagenome assembly using repeat graphs. **_Nat Methods_** 17, 1103–1110 **(2020)**. https://doi.org/10.1038/s41592-020-00971-x

- [`Canu`](https://github.com/marbl/canu) - (C++/C, v2.3, 2024.12)
    - Koren S, Walenz BP, Berlin K, et al. Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation. **_Genome Res_** 27, 722-736 **(2017)**. https://genome.cshlp.org/content/27/5/722

- [`FALCON`](https://github.com/PacificBiosciences/FALCON) - (Python, v0.3.0, 2018.8)

- [`miniasm`](https://github.com/lh3/miniasm) - (TeX/C, v0.3, 2018.7)
    - Li H. Minimap and miniasm: fast mapping and de novo assembly for noisy long sequences. **_Bioinformatics_** 32(14), 2103–2110 **(2016)**. https://doi.org/10.1093/bioinformatics/btw152

- [`wtdbg2`](https://github.com/ruanjue/wtdbg2) - (C, v2.5, 2019.09)
    - Ruan J, Li H. Fast and accurate long-read assembly with wtdbg2. **_Nat Methods_** 17, 155–158 **(2020)**. https://doi.org/10.1038/s41592-019-0669-3

- [`Trycycler`](https://github.com/rrwick/Trycycler) - (Python, v0.5.5, 2024.3)
    - Wick RR, Judd LM, Cerdeira LT, et al. Trycycler: consensus long-read assemblies for bacterial genomes. **_Genome Biol_** 22, 266 **(2021)**. https://doi.org/10.1186/s13059-021-02483-z

- [`NextDenovo`](https://github.com/Nextomics/NextDenovo) - (C, 2.5.2, 2023.3)
    - Hu J, Wang Z, Sun Z, et al. NextDenovo: an efficient error correction and accurate assembly tool for noisy long reads. **_Genome Biol_** 25, 107 **(2024)**. https://doi.org/10.1186/s13059-024-03252-4

- [`Autocycler`](https://github.com/rrwick/Autocycler) - (Rust, v0.4.0, 2025.4) - A tool for generating consensus long-read assemblies for bacterial genomes
    - Wick RR, Howden BP, Stinear TP. Autocycler: long-read consensus assembly for bacterial genomes. **bioRxiv** **(2025)**. https://doi.org/10.1101/2025.05.12.653612

- [`nanoMDBG`](https://github.com/GaetanBenoitDev/metaMDBG) - (C++, v1.1, 2024.12)
    - Benoit G, James R, Raguideau S, et al. High-quality metagenome assembly from nanopore reads with nanoMDBG. **bioRxiv** **(2025)**. https://doi.org/10.1101/2025.04.22.649928

- [`myloasm`](https://github.com/bluenote-1577/myloasm/) - (Rust, v0.1.0, 2025.5) - A new high-resolution long-read metagenome assembler for even noisy reads, documentaion: https://myloasm-docs.github.io

### Hybrid assembly
- [`idba/IDBA-Hybrid`](https://github.com/loneknightpy/idba) - (C++, v1.1.3, 2016.7)
    - Peng Y, Leung HCM, Yiu SM, Chin FYL. IDBA-UD: a de novo assembler for single-cell and metagenomic sequencing data with highly uneven depth. **_Bioinformatics_** 28(11), 1420–1428 **(2012)**. https://doi.org/10.1093/bioinformatics/bts174

- [`SPAdes/metaSPAdes`](https://github.com/ablab/spades) - (C++, v4.0.0, 2024.6)
    - Nurk S, Meleshko D, Korobeynikov A and Pevzner PA. metaSPAdes: a new versatile metagenomic assembler. **_Genome Res_** 27, 824-834 **(2017)**. https://doi.org/10.1101/gr.213959.116

- [`Aviary`](https://github.com/rhysnewell/aviary) - (Python, v0.9.2, 2024.9) - A hybrid assembly and MAG recovery pipeline (and more!)

- [`OPERA-MS`](https://github.com/CSB5/OPERA-MS) - (C++/Python, v0.8.3, 2020.6)
    - Bertrand D, Shaw J, Kalathiyappan M, et al. Hybrid metagenomic assembly enables high-resolution analysis of resistance determinants and mobile elements in human microbiomes. **_Nat Biotechnol_** 37, 937–944 **(2019)**. https://doi.org/10.1038/s41587-019-0191-2

- [`Unicycler`](https://github.com/rrwick/Unicycler) - (C++, v0.5.1, 2024.7)
    - Wick RR, Judd LM, Gorrie CL, Holt KE. Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads. **_PLoS Comput Biol_** 13(6), e1005595 **(2017)**. https://doi.org/10.1371/journal.pcbi.1005595

- [`HyLight`](https://github.com/LuoGroup2023/HyLight) - (C++/Python, v1.0.0, 2024.9)
    - Kang X, Zhang W, Li Y, et al. HyLight: Strain aware assembly of low coverage metagenomes. **_Nat Commun_** 15, 8665 **(2024)**. https://doi.org/10.1038/s41467-024-52907-0

### PB HiFi read assembly
- [`hifiasm-meta`](https://github.com/xfengnefx/hifiasm-meta/) - (C++, v0.3.2, 2024.9)
    - Feng X, Li H. Evaluating and improving the representation of bacterial contents in long-read metagenome assemblies. **_Genome Biol_** 25, 92 **(2024)**. https://doi.org/10.1186/s13059-024-03234-6

- [`metaMDBG`](https://github.com/GaetanBenoitDev/metaMDBG) - (C++, v1.1, 2024.12)
    - Benoit G, Raguideau S, James R, et al. High-quality metagenome assembly from long accurate reads with metaMDBG. **_Nat Biotechnol_** 42, 1378–1383 **(2024)**. https://doi.org/10.1038/s41587-023-01983-6

### Polish
- [`FMLRC2`](https://github.com/HudsonAlpha/fmlrc2) - (Rust, v0.1.8, 2022.7)
    - Mak QXC, Wick RR, Holt JM, Wang JR. Polishing De Novo Nanopore Assemblies of Bacteria and Eukaryotes With FMLRC2. **_Mol Biol Evol_** 40(3), msad048 **(2023)**. https://doi.org/10.1093/molbev/msad048

- [`pilon`](https://github.com/broadinstitute/pilon) - (Scala, v1.24, 2021.1)
    - Walker BJ, Abeel T, Shea T, et al. Pilon: An Integrated Tool for Comprehensive Microbial Variant Detection and Genome Assembly Improvement. **_PLoS ONE_** 9(11), e112963 **(2014)**. https://doi.org/10.1371/journal.pone.0112963

- [`NextPolish`](https://github.com/Nextomics/NextPolish) - (Rust, v1.4.1, 2022.7)
    - Hu J, Fan J, Sun Z, Liu S. NextPolish: a fast and efficient genome polishing tool for long-read assembly. **_Bioinformatics_** 36(7), 2253–2255 **(2020)**. https://doi.org/10.1093/bioinformatics/btz891

- [`Polypolish`](https://github.com/rrwick/Polypolish) - (Rust, v0.6.0, 2024.1)
    - Bouras G, Judd L, Edwards R, et al. How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies. **_Microbial Genomics_** 10, 6 **(2024)**. https://doi.org/10.1099/mgen.0.001254

- [`Pypolca`](https://github.com/gbouras13/pypolca) - (Python, v0.3.1, 2024.2)
    - Bouras G, Judd L, Edwards R, et al. How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies. **_Microbial Genomics_** 10, 6 **(2024)**. https://doi.org/10.1099/mgen.0.001254

- [`MaSuRCA/POLCA`](https://github.com/alekseyzimin/masurca) - (M4/Perl, v4.1.1, 2024.2) 
    - Zimin AV, Salzberg SL. The genome polishing tool POLCA makes fast and accurate corrections in genome assemblies. **_PLoS Comput Biol_** 16(6), e1007981 **(2020)**. https://doi.org/10.1371/journal.pcbi.1007981

- [`Racon`](https://github.com/lbcb-sci/racon) - (C++, v1.5.0, 2021.12)
    - Vaser R, Sović I, Nagarajan N and Šikić M. Ultrafast consensus module for raw de novo genome assembly of long uncorrected reads. **_Genome Res_** 27, 737-746 **(2017)**. http://www.genome.org/cgi/doi/10.1101/gr.214270.116

- [medaka]()
- [dorado]() - 

- [`DeepPolisher`](https://github.com/google/deeppolisher)
    - Mastoras M, Asri M, Brambrink L, et al. Highly accurate assembly polishing with DeepPolisher. **_Genome Res_** **(2025)**. https://doi.org/10.1101/gr.280149.124

### Genome size estimation
- [`LRGE`](https://github.com/mbhall88/lrge) - (Rust, v0.1.3, 2024.12) - Pure
    - Hall MB, Coin LJM. Genome size estimation from long read overlaps. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.11.27.625777

### Other elements assmebly
- [`Chopper`](https://github.com/kschimke/Chopper) - (Python) - NoRelease
    - Schimke KD, Vollmers C. Sequencing complete plasmids on Oxford Nanopore Technology Sequencers using R2C2 and Chopper. **bioRxiv** **(2025)**. https://doi.org/10.1101/2025.01.16.633418

- [`PlasMAAG`](https://github.com/RasmussenLab/vamb/tree/vamb_n2v_asy/workflow_PlasMAAG) - (snakemake) - part of VAMB
    - Líndez PP, Danielsen LS, Kovačić I, et al. Accurate plasmid reconstruction from metagenomics data using assembly-alignment graphs and contrastive learning. **bioRxiv** **(2025)**. https://doi.org/10.1101/2025.02.26.640269

### Assembly improvement
- [`COBRA`](https://github.com/linxingchen/cobra) - (Python, v1.2.3, 2024.2)
    - Chen L, Banfield JF. COBRA improves the completeness and contiguity of viral genomes assembled from metagenomes. **_Nat Microbiol_** 9, 737–750 **(2024)**. https://doi.org/10.1038/s41564-023-01598-2

## Binning
### Tools and workflows
- [`COCACOLA`](https://github.com/younglululu/COCACOLA) - (MATLAB/C++, 2015)
    - Lu YY, Chen T, Fuhrman JA, Sun F. COCACOLA: binning metagenomic contigs using sequence COmposition, read CoverAge, CO-alignment and paired-end read LinkAge, **_Bioinformatics_** 33(6), 791–798 **(2017)**. https://doi.org/10.1093/bioinformatics/btw290

- [`MetaBAT2`](https://bitbucket.org/berkeleylab/metabat) - (C++, v2.17, 2023.09)
    - Kang DD, Li F, Kirton E, et al. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. **_PeerJ_** 7, e7359 **(2019)**. https://doi.org/10.7717/peerj.7359

- [`MaxBin2`](https://sourceforge.net/projects/maxbin2/) - (C++, 2.2.7, 2020.06)
    - Wu YW, Simmons BA, Singer SW. MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. **_Bioinformatics_** 32(4), 605–607 **(2016)**. https://doi.org/10.1093/bioinformatics/btv638

- [`CONCOCT`](https://github.com/BinPro/CONCOCT) - (Python, v1.1.0, 2019.08)
    - Alneberg J, Bjarnason B, de Bruijn I, et al. Binning metagenomic contigs by coverage and composition. **_Nat Methods_** 11, 1144–1146 **(2014)**. https://doi.org/10.1038/nmeth.3103

- [`VAMB`](https://github.com/RasmussenLab/vamb) - (Python, v4.1.3, 2023.06)
    - Nissen JN, Johansen J, Allesøe RL, et al. Improved metagenome binning and assembly using deep variational autoencoders. **_Nat Biotechnol_** 39, 555–560 **(2021)**. https://doi.org/10.1038/s41587-020-00777-4
    
- [`TaxVAMB`](https://github.com/RasmussenLab/vamb) - (Python, v4.1.3, 2023.6)
    - Kutuzova S, Piera P, Nielsen KN, et al. Binning meets taxonomy: TaxVAMB improves metagenome binning using bi-modal variational autoencoder. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.10.25.620172

- [`UniteM`](https://github.com/dparks1134/UniteM) - (Python, 1.2.4, 2022.10)

- [`SemiBin`](https://github.com/BigDataBiology/SemiBin) - (Python, v1.5.1, 2023.03)
    - Pan S, Zhu C, Zhao XM, et al. A deep siamese neural network improves metagenome-assembled genomes in microbiome datasets across different environments. **_Nat Commun_** 13, 2326 **(2022)**. https://doi.org/10.1038/s41467-022-29843-y

- [`SemiBin2`](https://github.com/BigDataBiology/SemiBin) - (Python, v2.1.0, 2024.03) - Benchmark: https://github.com/BigDataBiology/SemiBin2_benchmark
    - Pan S, Zhao XM, Coelho LP. SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing. **_Bioinformatics_** 39(Suppl_1), 21–29 **(2023)**. https://doi.org/10.1093/bioinformatics/btad209

- [`Smeta`](https://github.com/YuhaoZhangwow/SMeta) - (C++, 2024.09) - NoReleaseTag
    - Zhang Y, Cheng M, Ning K. SMeta. a binning tool using single-cell sequences to aid in reconstructing species from metagenome accurately. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.08.25.609542
    - 
- [`Bin Chiken`](https://github.com/AroneyS/binchicken) - (Python, v0.12.6, 2024.12)
    - Aroney STN, Newell RJP, Tyson GW, Woodcroft BJ. Bin Chicken: targeted metagenomic coassembly for the efficient recovery of novel genomes. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.11.24.625082

- [`MetaCoAG`](https://github.com/metagentools/MetaCoAG) - (Python, v1.2.2, 2024.9)
    - Mallawaarachchi V, Lin Y. (2022). MetaCoAG: Binning Metagenomic Contigs via Composition, Coverage and Assembly Graphs. **_Research in Computational Molecular Biology RECOMB_** 13278  **(2022)**. https://doi.org/10.1007/978-3-031-04749-7_5

- [`GraphBin`](https://github.com/metagentools/GraphBin) - (Python, v1.7.4, 2024.8)
    - Mallawaarachchi V, Wickramarachchi A, Lin Y. GraphBin: refined binning of metagenomic contigs using assembly graphs. **_Bioinformatics_** 36(11), 3307–3313 **(2020)**. https://doi.org/10.1093/bioinformatics/btaa180

- [`GraphBin2`](https://github.com/metagentools/GraphBin2) - (Python, v1.3.3, 2024.9)
    - Mallawaarachchi VG, Wickramarachchi AS & Lin Y. Improving metagenomic binning results with overlapped bins using assembly graphs. **_Algorithms Mol Biol_** 16, 3 **(2021)**. https://doi.org/10.1186/s13015-021-00185-6

- [`MetaBinner`](https://github.com/ziyewang/MetaBinner) - (Python/Perl, v1.4.4, 2022.9)
    - Wang Z, Huang P, You R, et al. MetaBinner: a high-performance and stand-alone ensemble binning method to recover individual genomes from complex microbial communities. **_Genome Biol_** 24, 1 **(2023)**. https://doi.org/10.1186/s13059-022-02832-6

### Strain-level resolve
- [`inStrain`](https://github.com/MrOlm/inStrain) - (Python, v1.9.0, 2024.5)
    - Olm MR, Crits-Christoph A, Bouma-Gregson K, et al. inStrain profiles population microdiversity from metagenomic data and sensitively detects shared microbial strains. **_Nat Biotechnol_** 39, 727–736 **(2021)**. https://doi.org/10.1038/s41587-020-00797-0

- [`Floria`](https://github.com/bluenote-1577/floria) - (Rust, v0.0.1, 2024.1)
    - Shaw J, Gounot JS, Chen H, et al. Floria: fast and accurate strain haplotyping in metagenomes. **_Bioinformatics_** 40(Suppl_1), 30–38 **(2024)**. https://doi.org/10.1093/bioinformatics/btae252

- [`Lorikeet`](https://github.com/rhysnewell/Lorikeet) - (Rust, v0.8.2, 2023.12) - Strain resolver for metagenomics 

- [`STRONG`](https://github.com/chrisquince/STRONG) - (Python, v0.0.1-beta, 2021.5)
    - Quince C, Nurk S, Raguideau S, et al. STRONG: metagenomics strain resolution on assembly graphs. **_Genome Biol_** 22, 214 **(2021)**. https://doi.org/10.1186/s13059-021-02419-7

- [`Strainberry`](https://github.com/rvicedomini/strainberry) - (Python, v1.1, 2021.5)
    - Vicedomini R, Quince C, Darling AE, et al. Strainberry: automated strain separation in low-complexity metagenomes using long reads. **_Nat Commun_** 12, 4485 **(2021)**. https://doi.org/10.1038/s41467-021-24515-9

- [`Strainy`](https://github.com/katerinakazantseva/strainy) - (Python, v1.1, 2024.9)
    - Kazantseva E, Donmez A, Frolova M, et al. Strainy: phasing and assembly of strain haplotypes from long-read metagenome sequencing. **_Nat Methods_** 21, 2034–2043 **(2024)**. https://doi.org/10.1038/s41592-024-02424-1

- [`VStrains`](https://github.com/metagentools/VStrains) - (Python, v1.1.0, 2023.2)
    - Luo R, Lin Y. VStrains: De Novo Reconstruction of Viral Strains via Iterative Path Extraction From Assembly Graphs. **bioRxiv** **(2023)**. https://doi.org/10.1101/2022.10.21.513181

- [`devider`](https://github.com/bluenote-1577/devider) - (Rust, v0.0.1, 2024.10)
    - Shaw J, Boucher C, Yu YW, et al. devider: long-read reconstruction of many diverse haplotypes. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.11.05.621838

- [`ChronoStrain`](https://github.com/gibsonlab/chronostrain) - (Python, v0.6.0, 2025.4)
    - Kim Y, Worby CJ, Acharya S, et al. Longitudinal profiling of low-abundance strains in microbiomes with ChronoStrain. **_Nat Microbiol_** 10, 1184–1197 **(2025)**. https://doi.org/10.1038/s41564-025-01983-z

### MAG improvement
- [`Circlator`](https://github.com/sanger-pathogens/circlator) - (Python, v1.5.5, 2020.10)
    - Hunt M, Silva ND, Otto TD, et al. Circlator: automated circularization of genome assemblies using long sequencing reads. **_Genome Biol_** 16, 294 **(2015)**. https://doi.org/10.1186/s13059-015-0849-0

- [`gapFinisher`](https://github.com/kammoji/gapFinisher) - 2019
    - Kammonen JI, Smolander O-P, Paulin L, et al. gapFinisher: A reliable gap filling pipeline for SSPACE-LongRead scaffolder output. **_PLoS ONE_** 14(9), e0216885 **(2019)**. https://doi.org/10.1371/journal.pone.0216885

- [`FGAP`](https://github.com/pirovc/fgap) - (MATLAB, v1.8.1, 2017.11)
    - Piro VC, Faoro H, Weiss VA, et al. FGAP: an automated gap closing tool. **_BMC Res Notes_** 7, 371 **(2014)**. https://doi.org/10.1186/1756-0500-7-371

- [`MaSuRCA/samba.sh`](https://github.com/alekseyzimin/masurca)
    - Zimin AV, Salzberg SL. The SAMBA tool uses long reads to improve the contiguity of genome assemblies. **_PLoS Comput Biol_** 18(2), e1009860 **(2022)**. https://doi.org/10.1371/journal.pcbi.1009860

### MAG assessment
- [`BUSCO`](https://busco.ezlab.org) - (Python, v5.8.1, 2024.10) 
    - Manni M, Berkeley MR, Seppey M, et al. BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. **_Mol Biol Evol_** 38(10), 4647–4654 **(2021)**. https://doi.org/10.1093/molbev/msab199

- [`CheckM`](https://github.com/Ecogenomics/CheckM) - (Python, v1.2.3, 2024.06)
    - Parks DH, Imelfort M, Skennerton C, et al. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. **_Genome Res_** 25, 1043-1055 **(2015)**. https://doi.org/10.1101/gr.186072.114

- [`CheckM2`](https://github.com/chklovski/CheckM2) - (Python, v1.0.2, 2023.05)
    - Chklovski A, Parks DH, Woodcroft BJ et al. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. **_Nat Methods_** 20, 1203–1212 **(2023)**. https://doi.org/10.1038/s41592-023-01940-w

- [`RefineM`](https://github.com/dparks1134/RefineM) - (Python, v0.1.2, 2020.11) - Unsupported
    - Parks DH, Rinke C, Chuvochina M, et al. Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life. **_Nat Microbiol_** 2, 1533–1542 **(2017)**. https://doi.org/10.1038/s41564-017-0012-7

- [`MAGpurify`](https://github.com/snayfach/MAGpurify) - (Python, v2.1.2, 2020.03)
    - Nayfach S, Shi ZJ, Seshadri R, et al. New insights from uncultivated genomes of the global human gut microbiome. **_Nature_** 568, 505–510 **(2019)**. https://doi.org/10.1038/s41586-019-1058-x

- [`GUNC`](https://github.com/grp-bork/gunc) - (Python, v1.0.6, 2023.11)
    - Orakov A, Fullam A, Coelho LP, et al. GUNC: detection of chimerism and contamination in prokaryotic genomes. **_Genome Biol_** 22, 178 **(2021)**. https://doi.org/10.1186/s13059-021-02393-0

- [`DFAST_QC`](https://github.com/nigyta/dfast_qc) - (Python, v1.0.5, 2024.09)
    - Elmanzalawi M, Fujisawa T, Mori H, et al. DFAST_QC: quality assessment and taxonomic identification tool for prokaryotic Genomes. **_BMC Bioinformatics_** 26, 3 **(2025)**. https://doi.org/10.1186/s12859-024-06030-y

- [`dRep`](https://github.com/MrOlm/drep) - (Python, v3.4.2, 2023.02)
    - Olm MR, Brown CT, Brooks B, Banfield JF. dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes through de-replication. **_The ISME Journal_** 11(12), 2864–2868 **(2017)**. https://doi.org/10.1038/ismej.2017.126

- [`galah`](https://github.com/wwood/galah) - (Rust, v0.4.2, 2024.09) - More scalable dereplication for metagenome assembled genomes

- ['skDER'](https://github.com/raufs/skDER) - (Python, v1.3.3, 2025.7)
    - Salamzade R, Kottapalli A, Kalan LR. skDER and CiDDER: two scalable approaches for microbial genome dereplication. **_Microbial Genomics_** 11(7) **(2025)**. https://doi.org/10.1099/mgen.0.001438

- [`DAS Tool`](https://github.com/cmks/DAS_Tool) - (R/Ruby, v1.1.7, 2024.01)
    - Sieber CMK, Probst AJ, Sharrar A, et al. Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. **_Nat Microbiol_** 3, 836–843 **(2018)**. https://doi.org/10.1038/s41564-018-0171-1

- [`Binette`](https://github.com/genotoul-bioinfo/Binette) - (Python, v1.0.3, 2024.9) - A fast and accurate binning refinement tool to constructs high quality MAGs from the output of multiple binning tools

- [`MAGqual`](https://github.com/ac1513/MAGqual) - (Python, v0.3.0, 2024.8)
    - Cansdale A, Chong JPJ. MAGqual: a stand-alone pipeline to assess the quality of metagenome-assembled genomes. **_Microbiome_** 12, 226 **(2024)**. https://doi.org/10.1186/s40168-024-01949-z

### Other types of genome recovery
- [`BEREN`](https://gitlab.com/benminch1/BEREN) - (Python) -  No Release Tag
    - Minch B, Moniruzzaman M. BEREN: A bioinformatic tool for recovering Giant viruses, Polinton-like Viruses, and Virophages in metagenomic data. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.10.09.617401

## Gene and element prediction
### ORF
- [`Prodigal`](https://github.com/hyattpd/Prodigal) - (C, v2.6.3, 2016.02)
    - Hyatt D, Chen GL, LoCascio PF, et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. **_BMC Bioinformatics_** 11, 119 **(2010)**. https://doi.org/10.1186/1471-2105-11-119

- [`Pyrodigal`](https://github.com/althonos/pyrodigal) - (Cython/Python, v3.5.2, 2024.9)
    - Larralde M. Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. **_J Open Source Software_** 7(72), 4296 **(2022)**. https://doi.org/10.21105/joss.04296

- [`Prokka`](https://github.com/tseemann/prokka) - (Perl, v1.14.5, 2019.11)
    - Seemann T. Prokka: rapid prokaryotic genome annotation. **_Bioinformatics_** 30(14), 2068–2069 **(2014)**. https://doi.org/10.1093/bioinformatics/btu153

### Non-coding RNA
- [`Barrnap`](https://github.com/tseemann/barrnap) - (Perl/Shell, 0.9, 2018.4) - BAsic Rapid Ribosomal RNA Predictor. Barrnap predicts the location of ribosomal RNA genes in genomes
    
- [`RNAmmer`](https://services.healthtech.dtu.dk/services/RNAmmer-1.2/) - v1.2

- [`tRNAscan-SE 2.0`](https://github.com/UCSC-LoweLab/tRNAscan-SE) - (C/Perl, v2.0.12, 2022.11)
    - Chan PP, Lin BY, Mak AJ, et al. tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. **_Nucleic Acids Res_** 49(16), 9077–9096 **(2021)**. https://doi.org/10.1093/nar/gkab688

- [`tRNAscan-SE`](https://trna.ucsc.edu/tRNAscan-SE/) - Web: https://trna.ucsc.edu/tRNAscan-SE
    - Chan PP, Lowe TM. tRNAscan-SE: Searching for tRNA Genes in Genomic Sequences. In: Kollmar, M. (eds) Gene Prediction. **_Methods in Molecular Biology_** 1962. Humana, New York, NY (2019). https://doi.org/10.1007/978-1-4939-9173-0_1

- [`PILER-CR`](http://www.drive5.com/pilercr/) - v1.06
    - Edgar RC. PILER-CR: Fast and accurate identification of CRISPR repeats. BMC Bioinformatics 8, 18 (2007). https://doi.org/10.1186/1471-2105-8-18

- [`pybarrnap`](https://github.com/moshi4/pybarrnap) - (Python, v0.5.0, 2024.3) - pybarrnap: Python implementation of barrnap [Computer software]

### Signal peptides
- [`SignalP`](https://github.com/fteufel/signalp-6.0) - (Python, v6, 2022.1) - online prediction: https://services.healthtech.dtu.dk/services/SignalP-6.0/
    - Teufel F, Almagro Armenteros JJ, Johansen AR, et al. SignalP 6.0 predicts all five types of signal peptides using protein language models. **_Nat Biotechnol_** 40, 1023–1025 **(2022)**. https://doi.org/10.1038/s41587-021-01156-3

- [`DeepLocPro`](https://github.com/Jaimomar99/deeplocpro) - web api: https://ku.biolib.com/deeplocpro
    - Moreno J, Nielsen H, Winther O, et al. Predicting the subcellular location of prokaryotic proteins with DeepLocPro. **_Bioinformatics_** 40(12), btae677 **(2024)**. https://doi.org/10.1093/bioinformatics/btae677

### Plasmid databases
- [`PIPdb`](https://nmdc.cn/pipdb) - plasmids in pathogens, web: https://nmdc.cn/pipdb, 2024.9 
    - Zhu Q, Chen Q, Gao S, et al. PIPdb: a comprehensive plasmid sequence resource for tracking the horizontal transfer of pathogenic factors and antimicrobial resistance genes. **_Nucleic Acids Res_** 53(D1), D169–D178 **(2025)**. https://doi.org/10.1093/nar/gkae952

- [`PlasmidScope`](https://plasmid.deepomics.org/) - Website: https://plasmid.deepomics.org - download database
    - Li Y, Feng X, Chen X, et al. PlasmidScope: a comprehensive plasmid database with rich annotations and online analytical tools. **_Nucleic Acids Res_** 53(D1), D179–D188 **(2025)**. https://doi.org/10.1093/nar/gkae930

## Taxonomy profile
### Profile read
- [`MetaPhlAn`](https://github.com/biobakery/MetaPhlAn) - (Python, v4.1.1, 2024.5)
    - Blanco-Míguez A, Beghini F, Cumbo F, et al. Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. **_Nat Biotechnol_** 41, 1633–1644 **(2023)**. https://doi.org/10.1038/s41587-023-01688-w

- [`StrainPhlAn`](https://github.com/biobakery/MetaPhlAn) - (Python, v4.1.1, 2024.5)
    - Truong DT, Tett A, Pasolli E, et al. Microbial strain-level population structure and genetic diversity from metagenomes. **_Genome Res_** 27, 626-638 **(2017)**. https://dx.doi.org/10.1101/gr.216242.116

- [kraken2](https://github.com/DerrickWood/kraken2) - (C++, v2.1.3, 2023.6)
    - Wood DE, Lu J & Langmead B. Improved metagenomic analysis with Kraken 2. **_Genome Biol_** 20, 257 **(2019)**. https://doi.org/10.1186/s13059-019-1891-0

- [`Metabuli`](https://github.com/steineggerlab/Metabuli) - (C++, v1.0.7, 2024.9)
    - Kim J, Steinegger M. Metabuli: sensitive and specific metagenomic classification via joint analysis of amino acid and DNA. **_Nat Methods_** 21, 971–973 **(2024)**. https://doi.org/10.1038/s41592-024-02273-y

- [`Metabuli-App`](https://github.com/steineggerlab/Metabuli-App) - (js, v1.0.1, 2025.1) - compatible with macOS/Windows/Linux
    - Lee SJ, Kim J, Mirdita M, et al. Easy and interactive taxonomic profiling with Metabuli App. **bioRxiv** **(2025)**. https://doi.org/10.1101/2025.03.10.642298

- [`sylph`](https://github.com/bluenote-1577/sylph) - (Rust, v0.6.1, 2024.4)
    - Shaw J, Yu YW. Rapid species-level metagenome profiling and containment estimation with sylph. **_Nat Biotechnol_** **(2024)**. https://doi.org/10.1038/s41587-024-02412-y

- [`Kaiju`](https://github.com/bioinformatics-centre/kaiju) - (C, v1.10.1, 2024.3)
    - Menzel P, Ng K & Krogh A. Fast and sensitive taxonomic classification for metagenomics with Kaiju. **_Nat Commun_** 7, 11257 **(2016)**. https://doi.org/10.1038/ncomms11257

- [`mOTUs-db`](https://motus-db.org/) - web
    - Dmitrijeva M, Ruscheweyh HJ, Feer L, et al. The mOTUs online database provides web-accessible genomic context to taxonomic profiling of microbial communities. **_Nucleic Acids Res_** 53(D1), D797–D805 **(2025)**. https://doi.org/10.1093/nar/gkae1004

- [`mOTUs`](https://github.com/motu-tool/mOTUs) - (Python, v3.1.0, 2023.10)
    - Ruscheweyh HJ, Milanese A, Paoli L, et al. Cultivation-independent genomes greatly expand taxonomic-profiling capabilities of mOTUs across various environments. **_Microbiome_** 10, 212 **(2022)**. https://doi.org/10.1186/s40168-022-01410-z

- [`melon`](https://github.com/xinehc/melon) - (Python, v0.2.4, 2024.12)
    - Chen X, Yin X, Shi X, et al. Melon: metagenomic long-read-based taxonomic identification and quantification using marker genes. **_Genome Biol_** 25, 226 **(2024)**. https://doi.org/10.1186/s13059-024-03363-y

- [`kun-peng`](https://github.com/eric9n/Kun-peng) - (Rust, v0.7.4, 2024.9)
    - Chen Q, Zhang B, Peng C, et al. Kun-peng: an ultra-memory-efficient, fast, and accurate pan-domain taxonomic classifier for all. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.12.19.629356

- [`MNBC`](https://github.com/ComputationalPathogens/MNBC) - (Java, v1.2, 2025.1)
    - Lu R, Dumonceaux T, Anzar M, et al. MNBC: a multithreaded Minimizer-based Naïve Bayes Classifier for improved metagenomic sequence classification. **_Bioinformatics_** 40(10), btae601 **(2024)**. https://doi.org/10.1093/bioinformatics/btae601

### Profile contig
- [`MEGAN6`](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/megan6/)
    - Huson DH, Albrecht B, Bağcı C, et al. MEGAN-LR: new algorithms allow accurate binning and easy interactive exploration of metagenomic long reads and contigs. **_Biol Direct_** 13, 6 **(2018)**. https://doi.org/10.1186/s13062-018-0208-7
    
- [`MMseq2/easy-taxonomy`](https://github.com/soedinglab/mmseqs2) - (C/C++, v15-6f452, 2023.10)
    - Mirdita M, Steinegger M, Breitwieser F, et al. Fast and sensitive taxonomic assignment to metagenomic contigs. **_Bioinformatics_** 37(18), 3029–3031 **(2021)**. https://doi.org/10.1093/bioinformatics/btab184
    
- [`CDKAM`](https://github.com/SJTU-CGM/CDKAM) - (C++/Perl, v1.1, 2020.11)
    - Bui VK, Wei C. CDKAM: a taxonomic classification tool using discriminative k-mers and approximate matching strategies. **_BMC Bioinformatics_** 21, 468 **(2020)**. https://doi.org/10.1186/s12859-020-03777-y

- [`MetaMaps`](https://github.com/DiltheyLab/MetaMaps) - (Perl/C++)
    - Dilthey AT, Jain C, Koren S, et al. Strain-level metagenomic assignment and compositional estimation for long reads with MetaMaps. **_Nat Commun_** 10, 3066 **(2019)**. https://doi.org/10.1038/s41467-019-10934-2

- [`BugSeq`](https://app.bugseq.com/academic) - Web Online, 可试用
    - Fan J, Huang S & Chorlton SD. BugSeq: a highly accurate cloud platform for long-read metagenomic analyses. **_BMC Bioinformatics_** 22, 160 **(2021)**. https://doi.org/10.1186/s12859-021-04089-5

- [`Vamb/Taxometer`](https://github.com/RasmussenLab/vamb) - (Python, v4.1.3, 2023.6)
    - Kutuzova S, Nielsen M, Piera P, et al. Taxometer: Improving taxonomic classification of metagenomics contigs. **_Nat Commun_** 15, 8357 **(2024)**. https://doi.org/10.1038/s41467-024-52771-y

- [`geNomad`](https://github.com/apcamargo/genomad) - (Python, v1.8.0, 2024.4)
    - Camargo AP, Roux S, Schulz F, et al. Identification of mobile genetic elements with geNomad. **_Nat Biotechnol_** 42, 1303–1312 **(2024)**. https://doi.org/10.1038/s41587-023-01953-y

- [`VITAP`](https://github.com/DrKaiyangZheng/VITAP) - (Python, v1.7.1, 2025.2)
    - Zheng K, Sun J, Liang Y, et al. VITAP: a high precision tool for DNA and RNA viral classification based on meta-omic data. **_Nat Commun_** 16, 2226 **(2025)**. https://doi.org/10.1038/s41467-025-57500-7

### MAG taxonomy
- [`GTDB-Tk`](https://github.com/ecogenomics/gtdbtk) - (Python, v1.7.0, 2021.10) - A tools for genome classification based on GTDB. Allocating multiple CPUs can cause a dramatic increase in memory usage. Therefore, the program is typically run with a single CPU. It usually requires large memory, and small nodes (eg, with 120Gb RAM) are generally insufficient.
    - Chaumeil PA, Mussig AJ, Hugenholtz P, et al. GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. **_Bioinformatics_** 36(6), 1925–1927 **(2020)**. https://doi.org/10.1093/bioinformatics/btz848

- [`GTDB-Tk 2`](https://github.com/Ecogenomics/GTDBTk) - (Python, v2.4.1, 2025.04)
    - Chaumeil PA, Mussig AJ, Hugenholtz P, et al. GTDB-Tk v2: memory friendly classification with the genome taxonomy database. **_Bioinformatics_** 38(23), 5315–5316 **(2022)**. https://doi.org/10.1093/bioinformatics/btac672

- [`tronko`](https://github.com/lpipes/tronko) - C
    - Pipes L, Nielsen R. A rapid phylogeny-based method for accurate community profiling of large-scale metabarcoding datasets. **_eLife_** 13, e85794 **(2024)**. https://doi.org/10.7554/eLife.85794

- [`kMetaShot`](https://github.com/gdefazio/kMetaShot) - (Python,2024.9)
    - Defazio G, Tangaro MA, Pesole G, et al. kMetaShot: a fast and reliable taxonomy classifier for metagenome-assembled genomes. **_Briefings in Bioinformatics_** 26(1), bbae680 **(2025)**. https://doi.org/10.1093/bib/bbae680

## Annotation
### Tools for annotation
- [`eggNOG-mapper v2`](https://github.com/eggnogdb/eggnog-mapper) - (Python, v2.1.12, 2023.8)
    - Cantalapiedra CP, Hernández-Plaza A, Letunic I, et al. eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale. **_Mol Biol Evol_** 38(12), 5825–5829 **(2021)**. https://doi.org/10.1093/molbev/msab293

- [`KofamKOALA`](https://www.genome.jp/ftp/tools/kofam_scan/) - (Ruby, v1.3.0, 2020.05) - A web version 
    - Aramaki T, Blanc-Mathieu R, Endo H, et al. KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold. **_Bioinformatics_** 36(7), 2251–2252 **(2020)**. https://doi.org/10.1093/bioinformatics/btz859

- [`BlastKOALA`, `GhostKOALA`](https://www.kegg.jp/ghostkoala/) - Web
    - Kanehisa M, Sato Y, and Morishima K。 BlastKOALA and GhostKOALA: KEGG tools for functional characterization of genome and metagenome sequences. **_J Mol Biol_** 428, 726-731 **(2016)**. https://doi.org/10.1016/j.jmb.2015.11.006

- [`Macrel`](https://github.com/BigDataBiology/macrel) - (Python, v1.5.0, 2024.9)
    - Santos-Júnior CD, Pan S, Zhao X, Coelho LP. Macrel: antimicrobial peptide screening in genomes and metagenomes. **_PeerJ_** 8, e10555 **(2020)**. https://doi.org/10.7717/peerj.10555

- [`RGI`](https://github.com/arpcard/rgi) - (Python, 6.0.3, 2023.9) - This tool is designed for resistome prediction from protein or nucleotide data, with a complementary database required.
    - Alcock BP, Huynh W, Chalil R, et al. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. **_Nucleic Acids Res_** 51(D1), D690–D699 **(2023)**. https://doi.org/10.1093/nar/gkac920
    
- [`MetaCerberus`](https://github.com/raw-lab/metacerberus) - (Python, v1.4.0, 2024.8)
    - Figueroa JL III, Dhungel E, Bellanger M, et al. MetaCerberus: distributed highly parallelized HMM-based processing for robust functional annotation across the tree of life. **_Bioinformatics_** 40(3), btae119 **(2024)**. https://doi.org/10.1093/bioinformatics/btae119

- [`GMSC-mapper`](https://github.com/BigDataBiology/GMSC-mapper) - (Python, v0.1.0, 2024.04)
    - Duan Y, Santos-Júnior CD, Schmidt TS, et al. A catalog of small proteins from the global microbiome. **_Nat Commun_** 15, 7563 **(2024)**. https://doi.org/10.1038/s41467-024-51894-6

- [`Bakta`](https://github.com/oschwengers/bakta) - (Python, v1.10.1, 2024.11) - web: https://bakta.computational.bio
    - Schwengers O, Jelonek L, Dieckmann MA, et al. Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. **_Microbial Genomics_** 7, 11 **(2021)**. https://doi.org/10.1099/mgen.0.000685

- [`DRAM`](https://github.com/WrightonLabCSU/DRAM) - (Python, v1.5.0, 2024.1)
    - Shaffer M, Borton MA, McGivern BB, et al. DRAM for distilling microbial metabolism to automate the curation of microbiome function. **_Nucleic Acids Res_** 48(16), 8883–8900 **(2020)**. https://doi.org/10.1093/nar/gkaa621

- [`graftM`](https://github.com/geronimp/graftM) - (Python, v0.14.0, 2022.5)
    - Boyd JA, Woodcroft BJ, Tyson GW. GraftM: a tool for scalable, phylogenetically informed classification of genes within metagenomes. **_Nucleic Acids Res_** 46(10), e59 **(2018)**. https://doi.org/10.1093/nar/gky174

- [`CD-Search`](http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi) - Web
    - Marchler-Bauer A, Bryant SH. CD-Search: protein domain annotations on the fly. **_Nucleic Acids Res_** 32(suppl_2), W327–W331 **(2004)**. https://doi.org/10.1093/nar/gkh454

- [`HMMER3`](https://github.com/EddyRivasLab/hmmer) - (C, v3.4, 2023.8)
    - Eddy SR. Accelerated Profile HMM Searches. **_PLoS Comput Biol_** 7(10), e1002195 **(2011)**. https://doi.org/10.1371/journal.pcbi.1002195

- [`DeepARG`](https://github.com/gaarangoa/deeparg) - (Python2.7, 2023.11)
    - Arango-Argoty G, Garner E, Pruden A, et al. DeepARG: a deep learning approach for predicting antibiotic resistance genes from metagenomic data. **_Microbiome_** 6, 23 **(2018)**. https://doi.org/10.1186/s40168-018-0401-z

- [`ARGs-OAP`](https://github.com/xinehc/args_oap) - (Python, v3.2.4, 2023.10)
    - Yang Y, Jiang X, Chai B, et al. ARGs-OAP: online analysis pipeline for antibiotic resistance genes detection from metagenomic data using an integrated structured ARG-database. **_Bioinformatics_** 32(15), 2346–2351 **(2016)**. https://doi.org/10.1093/bioinformatics/btw136

- [`ARGs-OAP3.0`](https://github.com/xinehc/args_oap) - (Python, v3.2.4, 2023.10)
    - Yin X, Zheng X, Li L, et al. ARGs-OAP v3.0: Antibiotic-Resistance Gene Database Curation and Analysis Pipeline Optimization. **_Engineering_** 27, 234-241 **(2023)**. https://doi.org/10.1016/j.eng.2022.10.011

- [`argNorm`](https://github.com/BigDataBiology/argNorm) - (Python, v0.7.0, 2025.3)
    - Perovic SU, Ramji V, Chong H, et al. argNorm: normalization of antibiotic resistance gene annotations to the Antibiotic Resistance Ontology (ARO). **_Bioinformatics_** 41(5), btaf173 **(2025)**. https://doi.org/10.1093/bioinformatics/btaf173

- [`OrthoLoger`](https://orthologer.ezlab.org) - (bash, v3.5.0, 2024.10)
    - Kuznetsov D, Tegenfeldt F, Manni M, et al. OrthoDB v11: annotation of orthologs in the widest sampling of organismal diversity. **_Nucleic Acids Res_** 51(D1), D445–D451 **(2023)**. https://doi.org/10.1093/nar/gkac998

- [`nail`](https://github.com/TravisWheelerLab/nail) - (Rust, v0.2.0, 2024.7)
    - Roddy JW, Rich DH, Wheeler TJ. nail: software for high-speed, high-sensitivity protein sequence annotation. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.01.27.577580

- [`FastOMA`](https://github.com/DessimozLab/FastOMA/) - (Python, v0.3.4, 2024.9)
    - Majidian S, Nevers Y, Yazdizadeh Kharrazi A, et al. Orthology inference at scale with FastOMA. **_Nat Methods_** 22, 269–272 **(2025)**. https://doi.org/10.1038/s41592-024-02552-8

- [`BenchAMRking`](https://erasmusmc-bioinformatics.github.io/benchAMRking/) - web workflow AAMR detection
    - Strepis N, Dollee D, Vrins D, et al. BenchAMRking: a Galaxy-based platform for illustrating the major issues associated with current antimicrobial resistance (AMR) gene prediction workflows. **_BMC Genomics_** 26, 27 **(2025)**. https://doi.org/10.1186/s12864-024-11158-5

- [`OMArk`](https://github.com/DessimozLab/OMArk) - (Python, v0.3.0, 2023.10)
    - Nevers Y, Warwick Vesztrocy A, Rossier V, et al. Quality assessment of gene repertoire annotations with OMArk. **_Nat Biotechnol_** 43, 124–133 **(2025)**. https://doi.org/10.1038/s41587-024-02147-w

- [`mettannotator`](https://github.com/EBI-Metagenomics/mettannotator) - (Nextflow, 2024.12, v1.4.0)
    - Gurbich TA, Beracochea M, De Silva NH, et al. mettannotator: a comprehensive and scalable Nextflow annotation pipeline for prokaryotic assemblies. **_Bioinformatics_** 41(2), btaf037 **(2025)**. https://doi.org/10.1093/bioinformatics/btaf037

- [`pseudofinder`](https://github.com/filip-husnik/pseudofinder) - (Python, v1.1.0, 2022.3)
    - Syberg-Olsen MJ, Garber AI, Keeling PJ, et al. Pseudofinder: Detection of Pseudogenes in Prokaryotic Genomes. **_Mol Biol Evol_** 39(7), msac153 **(2022)**. https://doi.org/10.1093/molbev/msac153

- [`AntiDeffenseFinder`](https://defensefinder.mdmlab.fr/) - Web Service
    - Tesson F, Huiting E, Wei L, et al. Exploring the diversity of anti-defense systems across prokaryotes, phages and mobile genetic elements. **_Nucleic Acids Res_** 53(1), gkae1171 **(2025)**. https://doi.org/10.1093/nar/gkae1171

- [`Genomic conttext`](https://github.com/bio-ontology-research-group/Genomic_context) - Python - a step by step tutorial
    - Toibazar D, Kulmanov M, Hoehndorf R. Context-based protein function prediction in bacterial genomes. **bioRxiv** **(2024)**. Context-based protein function prediction in bacterial genomes

- [`OMAnnotator`](https://github.com/DessimozLab/OMAnnotator), Python
    - Bates S, Dessimoz C, Nevers Y. OMAnnotator: a novel approach to building an annotated consensus genome sequence. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.12.04.626846


- [`MHCScan`](https://github.com/Arkadiy-Garber/MHCscan) - Python
    - Garber AI, Nealson KH and Merino N. Large-scale prediction of outer-membrane multiheme cytochromes uncovers hidden diversity of electroactive bacteria and underlying pathways. **_Front Microbiol_** 15, 1448685 **(2024)**. https://doi.org/10.3389/fmicb.2024.1448685

- [`DRAMMA`](https://github.com/burstein-lab/DRAMMA) - Python
    - Rannon E, Shaashua S & Burstein D. DRAMMA: a multifaceted machine learning approach for novel antimicrobial resistance gene detection in metagenomic data. **_Microbiome_** 13, 67 **(2025)**. https://doi.org/10.1186/s40168-025-02055-4

- [`ChroQueTas`](https://github.com/nmquijada/ChroQueTas) - (Shell, v0.4.2, 2024.10) - Fungi, work with FungAMR
    - Bédard C, Pageau A, Fijarczyk A, et al. FungAMR: A comprehensive portrait of antimicrobial resistance mutations in fungi. **bioRxiv** **(2025)**. https://doi.org/10.1101/2024.10.07.617009

- [`MICROPHERRET`](https://github.com/BizzoTL/MICROPHERRET/) - Python
    - Bizzotto E, Fraulini S, Zampieri G, et al. MICROPHERRET: MICRObial PHEnotypic tRait ClassifieR using Machine lEarning Techniques. **_Environmental Microbiome_** 19, 58 **(2024)**. https://doi.org/10.1186/s40793-024-00600-6

### Databases for annotation
- [`eggNOG 6.0`](http://eggnog6.embl.de/) - (AA Genes, v6.0, 2022.09)
    - Hernández-Plaza A, Szklarczyk D, Botas J, et al. eggNOG 6.0: enabling comparative genomics across 12 535 organisms. **_Nucleic Acids Res_** 51(D1), D389–D394 **(2023)**. https://doi.org/10.1093/nar/gkac1022

- [`KEGG`](https://www.genome.jp/kegg/) - (AA Genes, 111.0, 2024.08)
    - Kanehisa M, Sato Y, Kawashima M, et al. KEGG as a reference resource for gene and protein annotation. **_Nucleic Acids Res_** 44(D1), D457–D462 **(2016)**. https://doi.org/10.1093/nar/gkv1070

- [`CAZy`](http://www.cazy.org) - AA, web
    - Drula E, Garron ML, Dogan S, et al. The carbohydrate-active enzyme database: functions and literature. **_Nucleic Acids Res_** 50(D1), D571–D577 **(2022)**. https://doi.org/10.1093/nar/gkab1045

- [`CARD`](https://card.mcmaster.ca) -（ARG, v3.3.0, 2024.08）
    - Alcock BP, Huynh W, Chalil R, et al. CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. **_Nucleic Acids Res_** 51(D1), D690–D699 **(2023)**. https://doi.org/10.1093/nar/gkac920

- [`UniProt`](https://www.uniprot.org/) - (AA genes, v2024_5) 
    - The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. **_Nucleic Acids Res_** 51(D1), D523–D531 **(2023)** , https://doi.org/10.1093/nar/gkac1052
    - The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. **_Nucleic Acids Res_** 53(D1), D609–D617 **(2025)**. https://doi.org/10.1093/nar/gkae1010

- [`ISOSDB`](https://github.com/joshuakirsch/pseudoR) - (IS Sequences, v3, 2023.12)
    - Kirsch JM, Hryckowian A, Duerkop B. A metagenomics pipeline reveals insertion sequence-driven evolution of the microbiota. **_Cell Host Microbe_** 32(5), 739-754.e4 **(2024)**. https://doi.org/10.1016/j.chom.2024.03.005

- [`AMPSphere`](https://ampsphere.big-data-biology.org/home) - (AMP, v2022-03, 2022.03)
    - Santos-Júnior CD, Torres MDT, Duan Y, et al. Discovery of antimicrobial peptides in the global microbiome with machine learning. **_Cell_** 187(14), 3761-3778.e16 **(2024)**. https://doi.org/10.1016/j.cell.2024.05.013

- [`DRAMP`](http://dramp.cpu-bioinfor.org) - (AMP, v4.0, 2024.09)
    - Shi G, Kang X, Dong F, et al. DRAMP 3.0: an enhanced comprehensive data repository of antimicrobial peptides. **_Nucleic Acids Res_** 50(D1), D488–D496 **(2022)**. https://doi.org/10.1093/nar/gkab651

- [`GMSC`](https://gmsc.big-data-biology.org) - (smORF, v1.0, 2024.08)
    - Duan Y, Santos-Júnior CD, Schmidt TS et al. A catalog of small proteins from the global microbiome. **_Nat Commun_** 15, 7563 **(2024)**. https://doi.org/10.1038/s41467-024-51894-6

- [`TCDB`](https://www.tcdb.org) - (transporters, 2024-09, 2024.09)
    - Saier MH, Reddy VS, Moreno-Hagelsieb G, et al. The Transporter Classification Database (TCDB): 2021 update. **_Nucleic Acids Res_** 49(D1), D461–D467 **(2021)**. https://doi.org/10.1093/nar/gkaa1004

- [`VFDB`](http://www.mgc.ac.cn/VFs/main.htm) - (virulence factors, 2024.9)
    - Liu B, Zheng D, Zhou S, et al. VFDB 2022: a general classification scheme for bacterial virulence factors. **_Nucleic Acids Res_** 50(D1), D912–D917 **(2022)**. https://doi.org/10.1093/nar/gkab1107

- [`DoriC`](https://tubic.org/doric/home) -  (oriC, v12.1)
    - Dong MJ, Luo H, Gao F. DoriC 12.0: an updated database of replication origins in both complete and draft prokaryotic genomes. **_Nucleic Acids Res_** 51(D1), D117–D120 **(2023)**. https://doi.org/10.1093/nar/gkac964

- [`TIGRAFMs`](http://tigrfams.jcvi.org/cgi-bin/index.cgi) - (protein sequence, v15.0, 2014.9)
    - Haft DH, Loftus BJ, Richardson DL, et al. TIGRFAMs: a protein family resource for the functional identification of proteins, **_Nucleic Acids Res_** 29(1), 41–43 **(2001)**. https://doi.org/10.1093/nar/29.1.41

- [`Pfam`](https://www.ebi.ac.uk/interpro/download/Pfam/) - (protein sequences, v37.0, 2024.5)
    - Mistry J, Chuguransky S, Williams L, et al. Pfam: The protein families database in 2021. **_Nucleic Acids Res_** 49(D1), D412–D419 **(2021)**. https://doi.org/10.1093/nar/gkaa913

- [`Rfam`](https://rfam.org) - (RNA families, v15, 2024.9)
    - Kalvari I, Nawrocki EP, Ontiveros-Palacios N, et al. Rfam 14: expanded coverage of metagenomic, viral and microRNA families. **_Nucleic Acids Res_** 49(D1), D192–D200 **(2021)**. https://doi.org/10.1093/nar/gkaa1047

- [`PHI`](https://phi5.phi-base.org) - (Pathogen Host Interactions, v5.0, 2024.3)
    - Urban M, Cuzick A, Seager J, et al. PHI-base in 2022: a multi-species phenotype database for Pathogen–Host Interactions. **_Nucleic Acids Res_** 50(D1), D837–D847 **(2022)**. https://doi.org/10.1093/nar/gkab1037

- [`oriTDB`](https://bioinfo-mml.sjtu.edu.cn/oriTDB2/)
    - Liu G, Li X, Guan J, et al. oriTDB: a database of the origin-of-transfer regions of bacterial mobile genetic elements. **_Nucleic Acids Res_** 53(D1), D163–D168 **(2025)**. https://doi.org/10.1093/nar/gkae869

- [`COG`](https://www.ncbi.nlm.nih.gov/research/COG) - FTP: https://ftp.ncbi.nlm.nih.gov/pub/COG/
    - Galperin MY, Alvarez RV, Karamycheva S, et al. COG database update 2024. **_Nucleic Acids Res_** 53(D1), D356–D363 **(2025)**. https://doi.org/10.1093/nar/gkae983

- [`OrthoDB`](https://www.orthodb.org) - v11.0
    - Kuznetsov D, Tegenfeldt F, Manni M, et al. OrthoDB v11: annotation of orthologs in the widest sampling of organismal diversity. **_Nucleic Acids Res_** 51(D1), D445–D451 **(2023)**. https://doi.org/10.1093/nar/gkac998

- [`CAZyme3D`](https://pro.unl.edu/CAZyme3D/) - web
    - Shanmugam NRS, Yin Y. CAZyme3D: a database of 3D structures for carbohydrate-active enzymes. **_J Mol Biol_** 437, 15 **(2025)**. https://doi.org/10.1016/j.jmb.2025.169001

- [`S9BactDB`](http://caps.ncbs.res.in/S9BactDB) - S9 proteases, web
    - Nayak S, Sowdhamini R. S9BactDB: A database for S9 family of proteases in bacterial genomes. **bioRxiv** **(2025)**. https://doi.org/10.1101/2025.01.01.631042

- [`FungAMR`](https://github.com/Landrylab/FungAMR) - (Fungi AMR) - work with ChroQueTas
    - Bédard C, Pageau A, Fijarczyk A, et al. FungAMR: A comprehensive portrait of antimicrobial resistance mutations in fungi. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.10.07.617009

### Gennome content analysis
- [`doubletrouble`](https://github.com/almeidasilvaf/doubletrouble) - (R, v1.6.0, 2025.2)
    - Almeida-Silva F, Van de Peer Y. doubletrouble: an R/Bioconductor package for the identification, classification, and analysis of gene and genome duplications. **_Bioinformatics_** 41(2), btaf043 **(2025)**. https://doi.org/10.1093/bioinformatics/btaf043

- [`gyoza`](https://github.com/durr1602/gyoza) - (snakemake, v11.0.0, 2025.2)
    - Durand R, Pageau A, Landry CR. gyōza: a Snakemake workflow for modular analysis of deep-mutational scanning data. **bioRxiv** **(2025)**. https://doi.org/10.1101/2025.02.19.639168

### Protein structure analysis
- [`Chai-1`](https://github.com/chaidiscovery/chai-lab) - (Python, v0.0.1, 2024.09) - blog: Introducing Chai-1: Decoding the molecular interactions of life, URL: https://www.chaidiscovery.com/blog/introducing-chai-1

- [`AlphaFold2`](https://github.com/google-deepmind/alphafold) - (Python, v2.3.2, 2023.0)
    - Jumper J, Evans R, Pritzel A, et al. Highly accurate protein structure prediction with AlphaFold. **_Nature_** 596, 583–589 **(2021)**. https://doi.org/10.1038/s41586-021-03819-2

- [`AlphaFold3`](https://www.alphafoldserver.com/) - Web Service
    - Abramson J, Adler J, Dunger J, et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. **_Nature_** 630, 493–500 **(2024)**. https://doi.org/10.1038/s41586-024-07487-w

- [`DNAproDB`](https://dnaprodb.usc.edu/) - Online database
    - Mitra R, Cohen AS, Sagendorf JM, et al. DNAproDB: an updated database for the automated and interactive analysis of protein–DNA complexes. **_Nucleic Acids Res_** 53(D1), D396–D402 **(2025)**. https://doi.org/10.1093/nar/gkae970

### Elements 
- [`DeepInverton`](https://github.com/HUST-NingKang-Lab/DeepInverton) - Python - NoReleaseTag
    - Wen J, Zhang H, Chu D, et al. Deep learning revealed the distribution and evolution patterns for invertible promoters across bacterial lineages. **_Nucleic Acids Res_** 52(21), 12817–12830 **(2024)**. https://doi.org/10.1093/nar/gkae966

### Bacteriophage
- [`SatelliteFinder`](https://hub.docker.com/r/gempasteur/satellite_finder) - Docker, Galaxy Service: https://galaxy.pasteur.fr/root?tool_id=toolshed.pasteur.fr/repos/fmareuil/satellitefinder/SatelliteFinder/0.9
    - Moura de Sousa JA, Fillol-Salom A, Penadés JR, et al. Identification and characterization of thousands of bacteriophage satellites across bacteria. **_Nucleic Acids Res_** 51(6), 2759–2777 **(2023)**. https://doi.org/10.1093/nar/gkad123

- [`Jaeger`](https://github.com/MGXlab/Jaeger) - (Python, v1.3.30-alpha, 2024.8)
    - Wijesekara Y, Wu LY, Beeloo R, et al. Jaeger: an accurate and fast deep-learning tool to detect bacteriophage sequences. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.09.24.612722

- [`PhageGE`](http://jason-zhao.shinyapps.io/PhageGE_Update/) - Web Platform
    - Zhao J, Han J, Lin YW, et al. PhageGE: an interactive web platform for exploratory analysis and visualization of bacteriophage genomes. **_GigaScience_** 13, giae074 **(2024)**. https://doi.org/10.1093/gigascience/giae074

- [Phynteny](https://github.com/susiegriggo/Phynteny_transformer) - (Jupyter, v0.1.2, 2025.06)
    - Grigson SR, Bouras G, Papudeshi B, et al. Synteny-aware functional annotation of bacteriophage genomes with Phynteny. **bioRxiv** **(2025)**. https://doi.org/10.1101/2025.07.28.667340

## Metabolic construction
### Tools for metabolic analysis
- [`antiSMASH`](https://github.com/antismash/antismash) - (Python, v7.1.0.1, 2023.11)
    - Blin K, Shaw S, Augustijn HE, et al. antiSMASH 7.0: new and improved predictions for detection, regulation, chemical structures and visualisation. **_Nucleic Acids Res_** 51(W1), W46–W50 **(2023)**. https://doi.org/10.1093/nar/gkad344

- [`MelonnPan`](https://github.com/biobakery/melonnpan) - R, dev
    - Mallick H, Franzosa EA, Mclver LJ, et al. Predictive metabolomic profiling of microbial communities using amplicon or metagenomic sequences. **_Nat Commun_** 10, 3136 **(2019)**. https://doi.org/10.1038/s41467-019-10927-1

- [`HUMAnN 2.0`](https://github.com/biobakery/humann/tree/2.9) - (Python, 最新2版v2.8.2, 2020.4)
    - Franzosa EA, McIver LJ, Rahnavard G, et al. Species-level functional profiling of metagenomes and metatranscriptomes. **_Nat Methods_** 15, 962–968 **(2018)**. https://doi.org/10.1038/s41592-018-0176-y

- [`HUMAnN 3.0`](https://github.com/biobakery/humann) - (Python, 最新3版v3.9, 2024.2)
    - Beghini F, McIver LJ, Blanco-Míguez A, et al. Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. **_eLife_** 10, e65088 **(2021)**. https://doi.org/10.7554/eLife.65088

- [`HUMAnN 4.0`](https://github.com/biobakery/humann) - (Python, 4.0.0.alpha.1-final, 22024.7) - active dev

- [`gapseq`](https://github.com/jotech/gapseq) - (R, v1.3.1, 2024.8)
    - Zimmermann J, Kaleta C & Waschina S. gapseq: informed prediction of bacterial metabolic pathways and reconstruction of accurate metabolic models.**_Genome Biol_** 22, 81 **(2021)**. https://doi.org/10.1186/s13059-021-02295-1

- [`./gapseq pan`](https://github.com/jotech/gapseq) - (R, v1.3.1, 2024.8) - The source code has been integrated into gapseq and is accessible via `./gapseq pan`
    - De Bernardini N, Zampieri G, Campanaro S, et al. pan-Draft: automated reconstruction of species-representative metabolic models from multiple genomes. **_Genome Biol_** 25, 280 **(2024)**. https://doi.org/10.1186/s13059-024-03425-1

- [`Bactabolize`](https://github.com/kelwyres/Bactabolize) - (Python, v1.0.4, 2024.12)
    - Vezina B, Watts SC, Hawkey J, et al. Bactabolize is a tool for high-throughput generation of bacterial strain-specific metabolic models. **_eLife_** 12, RP87406 **(2023)**. https://doi.org/10.7554/eLife.87406.3

### Metabolic databases
- [`mVOC 4.0`](http://bioinformatics.charite.de/mvoc)
    - Kemmler E, Lemfack MC, Goede A, et al. mVOC 4.0: a database of microbial volatiles. **_Nucleic Acids Res_** 53(D1), D1692–D1696 **(2025)**. https://doi.org/10.1093/nar/gkae961

- [`SMC`](https://smc.jgi.doe.gov) - Web
    - Udwary DW, Doering DT, Foster B, et al. The secondary metabolism collaboratory: a database and web discussion portal for secondary metabolite biosynthetic gene clusters. **_Nucleic Acids Res_** 53(D1), D717–D723 **(2025)**. https://doi.org/10.1093/nar/gkae1060

## Comparative genomics
### AAI and ANI
- [`ANI calculator`](http://enve-omics.ce.gatech.edu/ani/) - web

- [`AAI calculator`](http://enve-omics.ce.gatech.edu/aai/) - web

- [`fastANI`](https://github.com/ParBLiSS/FastANI) - (C++, v1.34, 2023.07)
    - Jain C, Rodriguez-R LM, Phillippy AM, et al. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. **_Nat Commun_** 9, 5114 **(2018)**. https://doi.org/10.1038/s41467-018-07641-9

- [`Pyskani`](https://github.com/althonos/pyskani) - (Python, v0.1.3, 2024.12) - [`PyOrthoANI`](https://github.com/althonos/pyorthoani) - (Python, v0.7.0, 2025.2) - [`PyFastANI`](https://github.com/althonos/pyfastani) - (Python, v0.6.1, 2025.1)
    - Larralde M, Zeller G, Carroll LM. PyOrthoANI, PyFastANI, and Pyskani: a suite of Python libraries for computation of average nucleotide identity. **bioRxiv**. **(2025)**. https://doi.org/10.1101/2025.02.13.638148

- [`EzAAI`](https://github.com/endixk/ezaai) - (Java, v1.2.3, 2024.02)
    - Kim D, Park S & Chun J. Introducing EzAAI: a pipeline for high throughput calculations of prokaryotic average amino acid identity. **_J Microbiol_** 59, 476–480 **(2021)**. https://doi.org/10.1007/s12275-021-1154-0

- [`SynTracker`](https://github.com/leylabmpi/SynTracker) - (Python, v1.3.1, 2024.09)
    - Enav H, Paz I & Ley RE. Strain tracking in complex microbiomes using synteny analysis reveals per-species modes of evolution. **_Nat Biotechnol_** 43, 773–783 **(2025)**. https://doi.org/10.1038/s41587-024-02276-2

- [`skani`](https://github.com/bluenote-1577/skani) - (Rust, v0.2.2, 2024.07)
    - Shaw J, Yu, YW. Fast and robust metagenomic sequence comparison through sparse chaining with skani. **_Nat Methods_** 20, 1661–1665 **(2023)**. https://doi.org/10.1038/s41592-023-02018-3

- [`ANIclustermap`](https://github.com/moshi4/ANIclustermap) - (Python, v1.4.0, 2024.7) 
    - Shimoyama Y. ANIclustermap: A tool for drawing ANI clustermap between all-vs-all microbial genomes. [Computer software] (2022). 

- [`EvANI`](https://github.com/sinamajidian/EvANI) - (Shell/Python, Norelease)
    - Sina Majidian, Stephen Hwang, Mohsen Zakeri, Ben Langmead. EvANI benchmarking workflow for evolutionary distance estimation. **bioRxiv**. **(2025)**. https://doi.org/10.1101/2025.02.23.639716

### view comparative map
- [`gggenes`](https://github.com/wilkox/gggenes) - Draw gene arrow maps in ggplot2

- [`LoVis4u`](https://github.com/art-egorov/lovis4u) - (Python, v0.0.11, 2024.10)
    - Egorov AA, Atkinson GC. LoVis4u: Locus Visualisation tool for comparative genomics. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.09.11.612399

- [`gggenomes`](https://github.com/thackl/gggenomes) - gggenomes: A Grammar of Graphics for Comparative Genomics

- [`plasmapR`](https://github.com/BradyAJohnston/plasmapR) - Creating plasmid maps inside ggplot

- [`geneviewer`](https://github.com/nvelden/geneviewer) - (R, CRAN release, 2025.1) - An R package designed for drawing gene arrow maps

### HGT
- [`WAFFLE`](https://github.com/biobakery/waafle) - (Python, v1.1.0, 2024.10)
    - Hsu TY, Nzabarushimana E, Wong D, et al. Profiling lateral gene transfer events in the human microbiome using WAAFLE. **_Nat Microbiol_** 10, 94–111 **(2025)**. https://doi.org/10.1038/s41564-024-01881-w

- [`SCARAP`](https://github.com/swittouck/scarap) - (Python, v1.0.0, 2024.11)
    - Wittouck S, Eilers T, van Noort V, et al. SCARAP: scalable cross-species comparative genomics of prokaryotes. **_Bioinformatics_** 41(1), btae735 **(2025)**. https://doi.org/10.1093/bioinformatics/btae735

### SV
- [`SVbyEye`](https://github.com/daewoooo/SVbyEye) - R, active dev
    - Porubsky D, Guitart X, Yoo DA, et al. SVbyEye: A visual tool to characterize structural variation among whole genome assemblies. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.09.11.612418

- [`CompareM`](https://github.com/dparks1134/CompareM) - (Python, v0.1.2, 2020.12) - Unsupported. Last version relsease on Dec 31 2020

- [`gcSV`](https://github.com/hitbc/gcSV) - C/C++, human
    - Li G, Liu Y, Liu B, et al. gcSV: a unified framework for comprehensive structural variant detection. **bioRxiv** **(2025)**. https://doi.org/10.1101/2025.02.10.637589

- [`minipileup`](https://github.com/lh3/minipileup) - (C, v1.1, 2025.2) - Minipileup is a simple pileup-based variant caller. It takes a reference FASTA and one or multiple alignment BAM as input, and outputs a multi-sample VCF along with allele counts

## Visualization
### View MSA
- [`MSAplot`](https://github.com/mourisl/MSAplot) - Python, no release tag. Plot multiple sequence alignment (MSA),with Jupyter notes

- [`pyMSAviz`](https://github.com/moshi4/pyMSAviz) - (Python, v0.5.0, 2024.9)
    - Shimoyama Y. pyMSAviz: MSA visualization python package for sequence analysis [Computer software] (2022)

- [`ggmsa`](https://github.com/YuLab-SMU/ggmsa) - (R, v1.0.2, 2021.8) - Visualizing publication-quality multiple sequence alignment using ggplot2,, CRAN移除了

- [`IGV`](https://igv.org) 
    - Thorvaldsdóttir H, Robinson JT, Mesirov JP. Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration. **_Briefings in Bioinformatics_** 14(2), 178–192 **(2013)**. https://doi.org/10.1093/bib/bbs017

- [`pymsaploter`](https://github.com/orangeSi/pymsaploter) - plot Multiple Sequence Alignment (MSA) of Clustal by python package, nuc only

### View genome
- [`pyGenomeViz`](https://github.com/moshi4/pyGenomeViz) - (Python, v1.4.1, 2024.9) - A genome visualization python package for comparative genomics

- [`pyCircos`](https://github.com/ponnhide/pyCircos) - (Python, v0.3.0, 2022.4)

- [`pyCirclize`](https://github.com/moshi4/pyCirclize) - (Python, v1.7.1, 2024.9) - Shimoyama Y. (2022). pyCirclize: Circular visualization in Python [Computer software]

- [`circlize`](https://cran.r-project.org/web/packages/circlize/) - (R, v0.4.16, 2024.2)
    - Gu Z, Gu L, Eils R, et al. circlize implements and enhances circular visualization in R. **_Bioinformatics_** 30(19), 2811–2812 **(2014)**. https://doi.org/10.1093/bioinformatics/btu393

### View assemblies
- [`Bandage`](https://github.com/rrwick/Bandage) - (C++, v0.9.0, 2022.1)
    - Wick RR, Schultz MB, Zobel J, et al. Bandage: interactive visualization of de novo genome assemblies. **_Bioinformatics_** 31(20), 3350–3352, **(2015)**. https://doi.org/10.1093/bioinformatics/btv383

## Phylogenetics
### Build a tree
- [`IQ-TREE`](http://www.iqtree.org) - (Java, v2.3.6, 2024.8)
    - Minh BQ, Schmidt HA, Chernomor O, et al. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. **_Mol Biol Evol_** 37(5), 1530–1534 **(2020)**. https://doi.org/10.1093/molbev/msaa015

- [`iqtree -m`](http://www.iqtree.org/#download) - (Java, v2.3.6, 2024.8)
    - Kalyaanamoorthy S, Minh B, Wong T, et al. ModelFinder: fast model selection for accurate phylogenetic estimates. **_Nat Methods_** 14, 587–589 **(2017)**. https://doi.org/10.1038/nmeth.4285

- [`fasttree2`](http://www.microbesonline.org/fasttree/) - (C++, v2.1.11)
    - Price MN, Dehal PS, Arkin AP. FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments. **_PLoS ONE_** 5(3), e9490 **(2010)**. https://doi.org/10.1371/journal.pone.0009490

- [`PhyloPhlAn`](https://github.com/biobakery/phylophlan) - (Python, v3.1.1, 2024.03)
    - Asnicar F, Thomas AM, Beghini F, et al. Precise phylogenetic analysis of microbial isolates and genomes from metagenomes using PhyloPhlAn 3.0. **_Nat Commun_** 11, 2500 **(2020)**. https://doi.org/10.1038/s41467-020-16366-7

- [`MrBayes`](https://github.com/NBISweden/MrBayes) - C
    - Huelsenbeck JP, Ronquist F. MRBAYES: Bayesian inference of phylogenetic trees. **_Bioinformatics_** 17(8), 754–755 **(2001)**. https://doi.org/10.1093/bioinformatics/17.8.754

- [`MrBayes`](https://github.com/NBISweden/MrBayes) - (C, v3)
    - Ronquist F, Huelsenbeck JP. MrBayes 3: Bayesian phylogenetic inference under mixed models. **_Bioinformatics_** 19(12), 1572–1574 **(2003)**. https://doi.org/10.1093/bioinformatics/btg180

- [`MrBayes`](https://github.com/NBISweden/MrBayes) - (C, v3.2.7, 2020.1)
    - Ronquist F, Teslenko M, van der Mark P, et al. MrBayes 3.2: Efficient Bayesian Phylogenetic Inference and Model Choice Across a Large Model Space. **_Systematic Biology_** 61(3), 539–542, **(2012)**. https://doi.org/10.1093/sysbio/sys029

- [`trimAI`](https://github.com/inab/trimal) - (C++, v1.5.0, 2024.7)
    - Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. **_Bioinformatics_** 25(15), 1972–1973 **(2009)**. https://doi.org/10.1093/bioinformatics/btp348

- [`RAxML-NG`](https://github.com/amkozlov/raxml-ng) - (C++, v1.2.2, 2024.5)
    - Kozlov AM, Darriba D, Flouri T, et al. RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. **_Bioinformatics_** 35(21), 4453–4455 **(2019)**. https://doi.org/10.1093/bioinformatics/btz305

- [`getphylo`](https://github.com/drboothtj/getphylo) - Python
    - Booth TJ, Shaw S, Cruz-Morales P, et al. getphylo: rapid and automatic generation of multi-locus phylogenetic trees. **_BMC Bioinformatics_** 26, 21 **(2025)**. https://doi.org/10.1186/s12859-025-06035-1

- [`Poplar`](https://github.com/sandialabs/poplar) - Python 
    - Elizabeth Koning, Raga Krishnakumar. Poplar: A Phylogenetics Pipeline. bioRxiv (2024). https://doi.org/10.1101/2024.11.11.623070
 
 - [`WitChi`](https://github.com/stephkoest/witchi) - Python
    - Köstlbacher S, Panagiotou K, Tamarit D, Ettema TJG. WitChi: Efficient Detection and Pruning of Compositional Bias in Phylogenomic Alignments Using Empirical Chi-Squared Testing. **bioRxiv** **(2025)**. https://doi.org/10.1101/2025.07.14.663642

### View tree
- [`iTol`](https://itol.embl.de) - Web
    - Letunic I, Bork P. Interactive Tree of Life (iTOL) v6: recent updates to the phylogenetic tree display and annotation tool. **_Nucleic Acids Res_** 52(W1), W78–W82 **(2024)**. https://doi.org/10.1093/nar/gkae268

- [`ARB`](http://www.arb-home.de) - (需要转发, arb-7.0, 2021.09)
    - Ludwig W, Strunk O, Westram R, et al. ARB: a software environment for sequence data. **_Nucleic Acids Res_** 32(4), 1363–1371 **(2004)**. https://doi.org/10.1093/nar/gkh293

- [`FigTree`](https://github.com/rambaut/figtree/) - (Java, v1.4.5-pre, 2024.02)

- [`GraPhlAn`](https://github.com/biobakery/graphlan) - (Python, 1.1.3, 2020.06)
    - Asnicar F, Weingart G, Tickle TL, et al. Compact graphical representation of phylogenetic data and metadata with GraPhlAn. **_PeerJ_** 3, e1029 **(2015)**. https://doi.org/10.7717/peerj.1029

- [`phyTreeViz`](https://github.com/moshi4/phyTreeViz) - (Python, v0.2.0, 2024.1) - Shimoyama Y. (2023). phyTreeViz: Simple phylogenetic tree visualization python package [Computer software]

## Microbial diversity analysis
### Abundance
- [`coverM`](https://github.com/wwood/CoverM) - (Rust, v0.7.0, 2024.1)
    - Aroney STN, Newell R, Nissen J, et al. CoverM: read alignment statistics for metagenomics. **_Bioinformatics_** 41(4), btaf147 **(2025)**. https://doi.org/10.1093/bioinformatics/btaf147

- [`Fairy`](https://github.com/bluenote-1577/fairy) - (Rust, v0.5.5, 2024.5)
    - Shaw J, Yu Y. Fairy: fast approximate coverage for multi-sample metagenomic binning. **_Microbiome_** 12, 151 **(2024)**. https://doi.org/10.1186/s40168-024-01861-6

- [`Tronko`](https://github.com/lpipes/tronko) - C, NoRelese
    - Pipes L, Nielsen R. A rapid phylogeny-based method for accurate community profiling of large-scale metabarcoding datasets. **_eLife_** 13, e85794 **(2024)**. https://doi.org/10.7554/eLife.85794

- [`MUSET`](https://github.com/CamilaDuitama/muset) - (C++, v0.5.1, 2024.12)
    - Vicedomini R, Andreace F, Dufresne Y, et al. MUSET: set of utilities for constructing abundance unitig matrices from sequencing data. **_Bioinformatics_** 41(3), btaf054 **(2025)**. https://doi.org/10.1093/bioinformatics/btaf054

- [`SingleM Microbial Fraction(SMF)`](https://github.com/EisenRa/2024_soil_dark_matter_reply) - Python 
    - Eisenhofer R, Alberdi A, Woodcroft BJ. Quantifying microbial DNA in metagenomes improves microbial trait estimation. **_ISME Communications_** 4(1), ycae111 **(2024)**. https://doi.org/10.1093/ismeco/ycae111

- [`Lyrebird`](https://wwood.github.io/singlem/Lyrebird) - Python - A tool developed by Samuel Aroney

- [`SingleM`](https://github.com/wwood/singlem) - (Python, v0.19.0, 2025.5)
    - Woodcroft BJ, Aroney STN, Zhao R, et al. Comprehensive taxonomic identification of microbial species in metagenomic data using SingleM and Sandpiper. **_Nat Biotechnol_** (2025). https://doi.org/10.1038/s41587-025-02738-1

### Diversity
- [`q2-kmerizer`](https://github.com/bokulich-lab/q2-kmerizer)
    - Bokulich NA. Integrating sequence composition information into microbial diversity analyses with k-mer frequency counting. **_mSystems_** 10 **(2025)**. https://doi.org/10.1128/msystems.01550-24

### Network
- [`SpeSpeNet`](https://tbb.bio.uu.nl/SpeSpeNet) - v1.0 - Web: https://tbb.bio.uu.nl/SpeSpeNet
    - van Eijnatten AL, van Zon L, Manousou E, et al. SpeSpeNet: an interactive and user-friendly tool to create and explore microbial correlation networks. **_ISME Communications_** ycaf036 **(2025)**. https://doi.org/10.1093/ismeco/ycaf036

- [`LUPINE`](https://github.com/SarithaKodikara/LUPINE) - (R, NoReleaseTag) - tutorial: https://mixomics.org/2024/05/lupine/
    - Kodikara S, Lê Cao KA. Microbial network inference for longitudinal microbiome studies with LUPINE. **_Microbiome_** 13, 64 **(2025)**. https://doi.org/10.1186/s40168-025-02041-w

### Interactions
- [`PHLAME`](https://github.com/quevan/phlame) - (Python, NoRelease Tag)
    - Qu EB, Baker JS, Markey L, et al. Intraspecies associations from strain-rich metagenome samples. **bioRxiv**. **(2025)** https://doi.org/10.1101/2025.02.07.636498 

- [`iPHoP`](https://bitbucket.org/srouxjgi/iphop/) - (Python, v1.3.3, 2023.11)
    - Roux S, Camargo AP, Coutinho FH, et al. iPHoP: An integrated machine learning framework to maximize host prediction for metagenome-derived viruses of archaea and bacteria. **_PLoS Biol_** 21(4), e3002083 **(2023)**. https://doi.org/10.1371/journal.pbio.3002083

- [`Anpan`](https://github.com/biobakery/anpan) - R
    - Ghazi AR, Thompson KN, Bhosle A, et al. Quantifying Metagenomic Strain Associations from Microbiomes with Anpan. **bioRxiv**. **(2025)**. https://doi.org/10.1101/2025.01.06.631550

## Metatranscriptomics
### RNA quantification
- [`sailfish`](https://github.com/kingsfordgroup/sailfish) - (C++, v0.10.0, 2016.4)
    - Patro R, Mount S & Kingsford C Sailfish enables alignment-free isoform quantification from RNA-seq reads using lightweight algorithms. **_Nat Biotechnol_** 32, 462–464 **(2014)**. https://doi.org/10.1038/nbt.2862

- [`SortMeRNA`](https://github.com/sortmerna/sortmerna) - (C++/C, v4.3.7, 2024.04)
    - Kopylova E, Noé L, Touzet H. SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. **_Bioinformatics_** 28(24), 3211–3217 **(2012)**. https://doi.org/10.1093/bioinformatics/bts611

- [`salmon`](https://github.com/COMBINE-lab/salmon) - (C++, v1.10.1, 2023.03)
    - Patro R, Duggal G, Love M, et al. Salmon provides fast and bias-aware quantification of transcript expression. **_Nat Methods_** 14, 417–419 **(2017)**. https://doi.org/10.1038/nmeth.4197

- [`RiboDetector`](https://github.com/hzi-bifo/RiboDetector) - (Python, v0.3.1, 2024.1)
    - Deng ZL, Münch PC, Mreches R, et al. Rapid and accurate identification of ribosomal RNA sequences via deep learning. **_Nucleic Acids Res_** 50(10), e60 **(2022)**. https://doi.org/10.1093/nar/gkac112

- [`kallisto`](https://github.com/pachterlab/kallisto) - (C/C++, v0.51.1, 2024.9)
    - Bray N, Pimentel H, Melsted P, et al. Near-optimal probabilistic RNA-seq quantification. **_Nat Biotechnol_** 34, 525–527 **(2016)**. https://doi.org/10.1038/nbt.3519

- [`NanoCount`](https://github.com/a-slide/NanoCount) - (Python, v1.1.0, 2022.12)
    - Gleeson J, Leger A, Prawer YDJ, et al. Accurate expression quantification from nanopore direct RNA sequencing with NanoCount. **_Nucleic Acids Res_** 50(4), e19 **(2022)**. https://doi.org/10.1093/nar/gkab1129

- [`step by step shell and python code`](https://github.com/RyoMameda/workflow)
    - Mameda R, Bono H. Data-driven workflow for comprehensive gene expression analysis in complex. **bioRxiv**. **(2025)**. https://doi.org/10.1101/2025.01.17.632662

### RNA assembly
- [`SPAdes/rnaSPAdes`](https://github.com/ablab/spades) - (C++, v4.0.0, 2024.6)
    - Bushmanova E, Antipov D, Lapidus A, et al. rnaSPAdes: a de novo transcriptome assembler and its application to RNA-Seq data. **_GigaScience_** 8(9), giz100 **(2019)**. https://doi.org/10.1093/gigascience/giz100

- [`idba/IDBA-tran`](https://github.com/loneknightpy/idba) - (C++, v1.1.3, 2016.7)
    - Peng Y, Leung HCM, Yiu SM, et al. IDBA-tran: a more robust de novo de Bruijn graph assembler for transcriptomes with uneven expression levels, **_Bioinformatics_** 29(13), 326–334 **(2013)**. https://doi.org/10.1093/bioinformatics/btt219

### RNA SV
- [`Clair3-RNA`](https://github.com/HKU-BAL/Clair3-RNA) - (Python, 2025.1, v0.2.1)
    - Zheng Z, Yu X, Chen L, et al. Clair3-RNA: A deep learning-based small variant caller for long-read RNA sequencing data. **bioRxiv** **(2024)**. https://doi.org/10.1101/2024.11.17.624050

### Stats
- [`GraphPCA`](https://github.com/YANG-ERA/GraphPCA) - (Python, GraphPCA-746ef16, 2024.8)
    - Yang J, Wang L, Liu L, et al. GraphPCA: a fast and interpretable dimension reduction algorithm for spatial transcriptomics data. **_Genome Biol_** 25, 287 **(2024)**. https://doi.org/10.1186/s13059-024-03429-x

## Surveillance
- [`VirDetector`](https://github.com/NLKaiser/VirDetector) - Nextflow, NoRelease
    - Kaiser NL, Groschup MH, Sadeghi B. VirDetector: a bioinformatic pipeline for virus surveillance using nanopore sequencing. **_Bioinformatics_** 41(2), btaf029 **(2025)**. https://doi.org/10.1093/bioinformatics/btaf029

- [`Castaanet`](https://github.com/MultipathogenGenomics/castanet) - (Python, v8.1, 2025.1)
    - Mayne R, Secret S, Geoghegan C, et al. Castanet: a pipeline for rapid analysis of targeted multi-pathogen genomic data. **_Bioinformatics_** 40(10), btae591 **(2024)**. https://doi.org/10.1093/bioinformatics/btae591

##  Modifications
- [`Robin`](https://github.com/looselab/robin) - (Python, Robin-d2f9b3, 2025.1) - A package to run real time analysis of nanopore methylation data
    - Deacon S, Cahyani I, Holmes N, et al. ROBIN: A unified nanopore-based assay integrating intraoperative methylome classification and next-day comprehensive profiling for ultra-rapid tumor diagnosis. **_Neuro-Oncology_** noaf103 **(2025)**. https://doi.org/10.1093/neuonc/noaf103

- [`uncalled4`](https://github.com/skovaka/uncalled4) - (Python, v4.1.0, 2024.8)
    - Kovaka S, Hook PW, Jenike KM, et al. Uncalled4 improves nanopore DNA and RNA modification detection via fast and accurate signal alignment. **_Nat Methods_** 22, 681–691 **(2025)**. https://doi.org/10.1038/s41592-025-02631-4

## Pangenome related
- [`PanGraph`](https://github.com/neherlab/pangraph) - (Julia/Python, v0.7.3, 2023.11)
    - Noll N​, Molari M​, Shaw LP and Neher RA. PanGraph: scalable bacterial pan-genome graph construction. **_Microbial Genomics_** 9, 6 **(2023)**. https://doi.org/10.1099/mgen.0.001034

- [`PanGenie`](https://github.com/eblerjana/pangenie) - (C++, v3.1.0, 2024.3)
    - Ebler J, Ebert P, Clarke WE et al. Pangenome-based genome inference allows efficient and accurate genotyping across a wide spectrum of variant classes. **_Nat Genet_** 54, 518–525 **(2022)**. https://doi.org/10.1038/s41588-022-01043-w

- [`Panaroo`](https://github.com/gtonkinhill/panaroo) - (Python, v1.5.0, 2024.4)
    - Tonkin-Hill G, MacAlasdair N, Ruis C, et al. Producing polished prokaryotic pangenomes with the Panaroo pipeline. **_Genome Biol_** 21, 180 **(2020)**. https://doi.org/10.1186/s13059-020-02090-4

- [`CELEBRIMBOR`](https://github.com/bacpop/CELEBRIMBOR) - (Python, v0.1.2, 2024.7)
    - Hellewell J, Horsfield ST, von Wachsmann J, et al. CELEBRIMBOR: core and accessory genes from metagenomes. **_Bioinformatics_** 40(9), btae542 **(2024)**. https://doi.org/10.1093/bioinformatics/btae542

- [`PangeBlocks`](https://github.com/AlgoLab/pangeblocks) - snakemake
    - Avila Cartes J, Bonizzoni P, Ciccolella S, et al. PangeBlocks: customized construction of pangenome graphs via maximal blocks. **_BMC Bioinformatics_** 25, 344 **(2024)**. https://doi.org/10.1186/s12859-024-05958-5

- [`PGGB`](https://github.com/pangenome/pggb) - (shell, v0.7.2, 2024.10)
    - Garrison E, Guarracino A, Heumos S, et al. Building pangenome graphs. **_Nat Methods_** 21, 2008–2012 **(2024)**. https://doi.org/10.1038/s41592-024-02430-3

- [`GrAnnoT`](https://forge.ird.fr/diade/dynadiv/grannot) - (Python, v1.0.3, 2024.8)
    - Marthe N, Zytnicki M, Sabot F. GrAnnoT, a tool for effecient and reliable annotation transfer through pangenome graph. **bioRxiv** **(2025)**. https://doi.org/10.1101/2025.02.26.640337

- [`Pangene`](https://github.com/lh3/pangene) - (C/JS, v1.1(r2311), 2024.2)
    - Li H, Marin M, Farhat MR. Exploring gene content with pangene graphs. **_Bioinformatics_** 40(7), btae456 **(2024)**. https://doi.org/10.1093/bioinformatics/btae456

- [`Mumemto`](https://github.com/vikshiv/mumemto) - (C++/Python, v1.3.1, 2025.6)
    - Shivakumar VS, Langmead B. Mumemto: efficient maximal matching across pangenomes. **_Genome Biol_** 26, 169 **(2025)**. https://doi.org/10.1186/s13059-025-03644-0

## Contact
This repository was created and maintained by [Jie Li](https://github.com/shley4).

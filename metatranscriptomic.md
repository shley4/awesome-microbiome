## metatranscriptomics
### meta RNA quantification
细胞中，不管是真核还是原核细胞中，rRNA通常占总体遗传物质的95%以上。在测序前，我们通常会通过一些实验操作将其从样本中去除，例如polyA方法，beads方法等。然而，由于RNA自身的特性，导致测序时常测到rRNA，但这部分的数据通常是会影响我们定量的结果的，所以我们必须在进行下游分析前去除掉这些rRNA的read。本集合收录了相关的工具。
- [[sailfish](https://github.com/kingsfordgroup/sailfish)] - [C++, v0.10.0, 2016.4] - [2014.4, _Nat Biotechnol_] - [[Sailfish enables alignment-free isoform quantification from RNA-seq reads using lightweight algorithms.](https://doi.org/10.1038/nbt.2862)]
- [[SortMeRNA](https://github.com/sortmerna/sortmerna)] - [C++/C, v4.3.7, 2024.04] - [2012.10, _Bioinformatics_] - [[SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data.](https://doi.org/10.1093/bioinformatics/bts611)]
- [[salmon](https://github.com/COMBINE-lab/salmon)] - [C++, v1.10.1, 2023.03] - [2017.03, _Nat Methods_] - [[Salmon provides fast and bias-aware quantification of transcript expression.](https://doi.org/10.1038/nmeth.4197)]
- [[RiboDetector](https://github.com/hzi-bifo/RiboDetector)] - [Python, v0.3.1, 2024.1] - [2022.6, _Nucleic Acids Res_] - [[Rapid and accurate identification of ribosomal RNA sequences via deep learning](https://doi.org/10.1093/nar/gkac112)]
- [[kallisto](https://github.com/pachterlab/kallisto)] - [C/C++, v0.51.1, 2024.9] - [2016.4, _Nat Biotechnol_] - [[Near optimal probabilistic RNA-seq quantification.](https://doi.org/10.1038/nbt.3519)]

## 组装 | assembly
- [[SPAdes/rnaSPAdes](https://github.com/ablab/spades)] -  [] - [2019.9, _Gigascience_] - [[rnaSPAdes: a de novo transcriptome assembler and its application to RNA-Seq data.](https://doi.org/10.1093/gigascience/giz100)]
- [[idba/IDBA-tran](https://github.com/loneknightpy/idba)] - [C++, v1.1.3, 2016.7] - [2013.7, _Bioinformatics_] - [[IDBA-tran: a more robust de novo de Bruijn graph assembler for transcriptomes with uneven expression levels.](https://doi.org/10.1093/bioinformatics/btt219)]

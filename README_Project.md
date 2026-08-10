# The impact of genetic background on eclosion timing across Drosophila timeless genotypes

\# Author: Cetline Sailer



This repository contains the data, analysis scripts, and computational workflows used for the Bachelor's thesis investigating the genetic basis of eclosion rhythm in *Drosophila melanogaster*.

The project focuses on the effect of the **tim01 mutation** and investigates whether its phenotypic effect on eclosion timing depends on the genomic background.

\---

## Repository structure

```text
.
├── code/
│   ├── Rscripts/ 
│   ├── Snakefile
│   ├── config.CS\_2026-07-24.yaml
│   ├── config\_cetline\_test.yaml
│   └── run\_genotype\_tim\_01\_CS\_2026-07-24.sh
│
├── conda/
│   └── envs/
│
└── data/
    ├── ID\_locus\_allele\_primers.txt
    ├── bed/                   -> flank/del bed files of timeless gene
    ├── fastq/                 -> DNA sequencing data of F2 
    ├── parental\_img\_analysis/ -> P0 generation recordings (tim01 and wildtype GP-N19)
    ├── recomb1\_img\_analysis/  -> F2 generation recordings (without control)
    ├── recomb2\_img\_analysis/  -> F2 generation recordings (with WT\_P0 control)
    ├── refs/                  -> timeless reference sequence
    └── test\_data/
```

\---

# 


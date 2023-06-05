## Getting started
- Install through pip: `pip install hifieval`
- Install through conda: `conda install hifieval`
- Install one error correction/assembly tool: [hifiasm](https://github.com/chhylp123/hifiasm) for example
- Install [minimap2](https://github.com/lh3/minimap2)
```
# get test data
wget https://zenodo.org/record/7799845/files/ecoli.reads.fastq?download=1  # simulated raw reads
wget https://zenodo.org/record/7799845/files/ecoli.ref.fasta?download=1  # reference genome

# get error corrected reads
hifiasm -o ecoli.asm.hifiasm --primary -t 10 --write-ec ecoli.reads.fastq 2> ecoli.asm.hifiasm.log

# get alignment paf files
minimap2 -t 8 -cx map-hifi --secondary=no --paf-no-hit --cs ecoli.ref.fasta ecoli.reads.fastq > ecoli.raw.paf
minimap2 -t 8 -cx map-hifi --secondary=no --paf-no-hit --cs ecoli.ref.fasta ecoli.asm.hifiasm.ec.fa > ecoli.hifiasm.paf

# get evaluation files
hifieval.py -o ecoli.hifiasm -r ecoli.raw.paf -c ecoli.hifiasm.paf
```

## Users' Guide

```hifieval [options] -r <raw.paf> -c <corrected.paf>```

Hifieval is a tool to evaluate long-read error correction mainly with PacBio High-Fidelity Reads (HiFi reads). Use command `hifieval` to see available options.

The input of this tool takes in two .paf files: one is raw reads aligned to reference genome; the other is corrected reads aligned to reference genome. PAF is a text format describing the approximate mapping positions between two set of sequences. 

The paf file will encodes difference of sequence alignments in the short form, indication substitution, insertion, and deletion. The metrics of error correction are:
- OC: (over-correction) The errors appeared in corrected reads but not in raw reads 
- UC: (under-correction) The errors in raw reads that are still in corrected reads
- CC: (correct-correction) The errors that are in raw reads but not corrected reads

### General usage
- Examples of Error Correction (EC) tools to output error corrected reads
    - hifiasm: ```hifiasm -o <prefix> --write-ec -t32 <read_files> 2> <prefix>.log```
    - LJA:     ```lja -o <output_dir> --reads <reads_file> [--reads <reads_file2> …]```
    - Verkko:  ```verkko -d <output_dir> --hifi <reads_files>```

- If the EC tool produce HPC corrected reads, use [seqtk](https://github.com/lh3/seqtk) to perform homopolymer-compression (HPC) on raw reads and the reference: 
```seqtk hpc <file>```
- [Minimap2](https://lh3.github.io/minimap2/minimap2.html) is used to generate the paf files using the command, the --cs tag is required: 
```./minimap2 -t8 -cx map-hifi --secondary=no --paf-no-hit --cs <ref_fasta_file> <read_files>  > <prefix>.paf```
    
### Advanced features
On top of FPR and TPR for the corrections, errors in homopolymer (HP) regions can be further incorporated if the assembly tool does not perform HPC on the raw reads during the error correction step using the command:

```hifieval [options] -h <reference_file> -r <raw.paf> -c <corrected.paf>```

HP regions of different lengths are identified, and UC/OC that fall within these regions is calculated. Here the error rate is calculated by ${\verb|#|HP_{x\,with\,error}}/{\verb|#|HP_{x}}$. for HP with length $x$. However, since most of the assembly tools use HPC reads during their error correction step, HP evaluation is optional. 

### Output overview
1. summary.tsv: the most detailed summary of EC performance for any downstream analysis
    - contains 12 columns: readName, raw_mapped_chr, raw_start, raw_end, raw_mq, corrected_mapped_chr, corrected_start, corrected_end, corrected_mq,  num_oc, num_uc, num_cc
3. rdlvl.eval.tsv
    - counts how many corrected reads have 1 oc/uc, 2 oc/uc, etc. for each chromosome and all chromosomes
5. metric.eval.tsv
    - overall metrics for each chromosome and all chromosomes
7. hp.ErrorRate.tsv
    - contains the error rates for each length of the homopolymers

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

```python hifieval.py -g <ref.fa>|<ref.fq> -r <raw.paf> -c <corr.paf>```

Hifieval is a tool to evaluate long-read error correction mainly with PacBio High-Fidelity Reads (HiFi reads). 

The input of this tool takes in two .paf files: one is raw reads aligned to reference genome; the other is corrected reads aligned to reference genome. PAF is a text format describing the approximate mapping positions between two set of sequences. 

[Minimap2](https://lh3.github.io/minimap2/minimap2.html) is used to generate the paf files using the command: 
```./minimap2 -t8 -cx map-hifi --secondary=no --paf-no-hit --cs <ref_fasta_file> <read_files>  > aln.paf```

With the --cs tag, the paf file will encodes difference sequences in the short form, indication substitution, insertion, and deletion. The metrics of error correction are:
- OC: (over-correction) The errors appeared in corrected reads but not in raw reads 
- UC: (under-correction) The errors in raw reads that are still in corrected reads
- CC: (correct-correction) The errors that are in raw reads but not corrected reads

### General usage
- Error Correction tools
Hifiasm:
```hifiasm -o prefix.asm --write-ec -t32 read_files 2> prefix.log```

LJA:

```lja [options] -o <output_dir> --reads <reads_file> [--reads <reads_file2> …]```

verkko:

```verkko -d <output_dir> --hifi <reads_files>```
    
### Advanced features
On top of FPR and TPR for the corrections, errors in homopolymer (HP) regions can be further incorporated if the assembly tool does not perform homopolymer-compression (HPC) on the raw reads during the error correction step. HP regions of different lengths are identified, and UC/OC that fall within these regions is calculated. Here the error rate is calculated by ${\#HP_{x\,with\,error}}/{\#HP_{x}}$.} for HP with length $x$. However, since most of the assembly tools use HPC reads during their error correction step, HP evaluation is optional. 


### Output overview

## Limitations

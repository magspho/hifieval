## Getting started
- Install through pip: `pip install hifieval`
- Install through conda: `conda install hifieval`

```
wget test data

```


## Users' Guide

```python hifieval.py -g <ref.fa>|<ref.fq> -r <raw.paf> -c <corr.paf>```

Hifieval is a tool to evaluate long-read error correction mainly with PacBio High-Fidelity Reads (HiFi reads). 

The genome and reads dataset come from [CHM13v1.1](https://github.com/marbl/CHM13). The fastq.gz files contain 30x PacBio HiFi reads. A smaller dataset from [*E. coli*](https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2/) is also used for testing and read simulation.

The input of this tool takes two .paf files: one is raw reads vs ref genome; the other is corrected reads vs ref genome. PAF is a text format describing the approximate mapping positions between two set of sequences. 

[Minimap2](https://lh3.github.io/minimap2/minimap2.html) is used to generate the paf files using the command: 
```./minimap2 -cx map-hifi -secondary=no --paf-no-hit --cs <ref_fasta_file> <read_files>  aln.paf```

With the -cs tag, the paf file will encodes difference sequences in the short form, indication substitution, insertion, and deletion. The metrics of error correction are:
- oc: (over-correction) The errors appeared in corrected reads (corr) but not in raw reads (raw)
- uc: (under-correction) The errors in raw that are still in corr
- cc: (correct-correction) The errors that are in raw but not corr

On top of FPR and TPR for the corrections, errors in homopolymer (HP) regions is further incorporated , since most of the tools use homopolymer-compressed (HPC) reads during sequence alignment. HP regions of different lengths are identified and uc/oc’s that fall within these regions. The HP error rate of length i = (# HP of len(i) with error) / (total # HP of len(i)). However, as most of the tools compress all HP’s, HP evaluation is optional. A non-HP evaluation using seqtk hpc is still under progress.
### General usage
- Error Correction tools
Hifiasm:

	./hifiasm/hifiasm -o prefix.asm --write-ec -t32 read_files 2>prefix.log

LJA:

	./LJA/bin/lja [options] -o <output_dir> --reads <reads_file> [--reads <reads_file2> …]

verkko:

	verkko -d <output_dir> --hifi <reads_files>
	
This canu-based tool is really new and performs really well. The 0-step is HiFi correction, it has a few subsets (k-mer counts, overlaps, correction) and it gets cleaned up as soon as the hifi.corrected.fasta file exists. Usually, when snakemake tries to go back to an earlier step it's because something either got manually changed in the folder or the input data got updated/touched.
    
### Advanced features

### Algorithm overview
### Getting help

## Developers' Guide
## Limitations

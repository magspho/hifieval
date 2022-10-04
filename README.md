# hifieval

## Getting started

how to download

## Users' Guide

### General usage
- Error Correction tools
HiCanu: 
	./canu -correct -p <prefix> -d <output_dir> genomeSize=4.8m -pacbio <read_files>
	This command can perform error correction without specifying pacbio-hifi
Canu selects multiple error rate thresholds and selects the most appropriate one based on how many reads end up without overlaps at each threshold. By default, all needed top-level tasks are performed (-pacbio and -nanopore are assumed to be raw and untrimmed while -pacbio-hifi are assumed to be corrected and trimmed). It is possible to run exactly one task by specifying your read characteristics and a step name. These options can be useful if you want to correct reads once and try many different assemblies. In the documentation, it says “HiCanu has support for PacBio HiFi data by compressing homopolymers, correcting isolated errors, and masking systematic errors”, but I have yet succeeded to generate corrected reads using the option –pacbio-hifi
mdBG:
	target/release/rust-mdbg  read_files utils/magic_simplify prefix
	mdBG uses a partial-order alignment approach, where we align reads from the same genomic region iteratively and come up with a consensus sequence in minimizer-space. Therefore, there won't be a bijection between reads and their corrected counterparts: There are more than one read that goes into POA (although you can upper bound this number with an input parameter), and we discard all reads that are used in the POA with the consensus. 

Hifiasm:
	./hifiasm/hifiasm -o prefix.asm --write-ec -t32 read_files 2>prefix.log

LJA:
	./LJA/bin/lja [options] -o <output_dir> --reads <reads_file> [--reads <reads_file2> …]

verkko:
	verkko -d <output_dir> --hifi <reads_files>
	This canu-based tool is really new and performs really well. The 0-step is HiFi correction, it has a few subsets (k-mer counts, overlaps, correction) and it gets cleaned up as soon as the hifi.corrected.fasta file exists. Usually, when snakemake tries to go back to an earlier step it's because something either got manually changed in the folder or the input data got updated/touched.

CONSENT:
	./CONSENT-correct --in <read_files> --out <out_file>-type PB
	This tool is no longer maintained by the developers, I kept running into pointer errors, but others on the Github issue posts seem to succeed once in a while.
    
### Advanced features
            ▪ Working with >65535 CIGAR operations
            ▪ The cs optional tag
            ▪ Working with the PAF format
### Algorithm overview
### Getting help

## Developers' Guide
## Limitations

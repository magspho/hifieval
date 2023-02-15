#!/usr/bin/env python

import pandas as pd
import matplotlib as plt
import numpy as np

tool_compare = ["hifiasm-22", "hifiasm",'lja','verkko','hifiasmhpc']
genome = "chm13"

# eval metric
metric_dict = dict()
for tool in tool_compare:
    filename = genome + '.' + tool + '.metric.eval.tsv'
    metric_dict[tool] = pd.read_csv(filename, sep='\t')
    
rates = pd.DataFrame(columns = ['FDR', 'FNR'])
for tool, metric in metric_dict.items():
    tmp = metric.loc[metric['chrName'] == 'all'][['FDR','FNR']]
    tmp.index = [tool]
    rates = pd.concat([rates,tmp])

rates.to_csv('chm13.tools.metric.eval.tsv', sep = '\t')
    
ax = rates.plot(kind='bar')
ax.figure.savefig('chm13.tools.metric.eval.pdf')

# readlevel eval
readinfo_dict = dict()
for tool in tool_compare:
    filename = genome + '.' + tool + '.rdlvl.eval.tsv'
    readinfo_dict[tool] = pd.read_csv(filename, sep='\t')

allchr_readinfo = pd.DataFrame()
for tool, readinfo in readinfo_dict.items():
    tmp = readinfo.loc[(readinfo['chr'] == 'all_uc') | (readinfo['chr'] == 'all_oc')]
    tmp = pd.DataFrame(tmp.sum()).tail(-1)
#     print(tool)
#     print(tmp)
    allchr_readinfo = pd.concat([allchr_readinfo, tmp], axis = 1)

allchr_readinfo.index = allchr_readinfo.index.astype(int)
allchr_readinfo = allchr_readinfo.sort_index().replace(np.nan, 0).head(500)
allchr_readinfo.columns = tool_compare
allchr_readinfo.reset_index(inplace=True)
allchr_readinfo.to_csv('chm13.tools.rdlvl.eval.tsv', sep = '\t')

ax1 = allchr_readinfo.plot(y = ['hifiasm-22', 'hifiasm'], kind="bar", figsize=(15,5), logy=True, xticks=np.arange(1,500,10))
ax2 = allchr_readinfo.plot(y = ['hifiasmhpc', 'lja', 'verkko'], kind="bar", figsize=(15,8), logy=True, xticks=np.arange(1,400,10))

ax1.figure.savefig('chm13.hifiasm.rdlvl.eval.pdf')
ax2.figure.savefig('chm13.tools.rdlvl.eval.pdf')


# def main(argv):
#     opts, args = getopt.getopt(argv[1:],"o:h:br:c:", 
#                                ["prefix=","hp=","specbed","raw=","corrected="])
    
#     if len(opts) < 2:
#         print("Usage: hifieval.py  [options]")
#         print("Options:")
#         print("  -o STR      Output File Prefix")
#         print("  -h STR      FASTA file with reference genome for evaluation in homopolymer region")
#         print("  -b STR      BED file with specified regions for evaluation")
#         print("  -r STR      PAF file aligned between raw reads and reference genome")
#         print("  -c STR      PAF file aligned between corrected reads and reference genome")
#         print("Minimap2 command for generating PAF file:")
#         print("  minimap2 -t32 -cx map-hifi --secondary=no --paf-no-hit --cs <reference genome file> <reads file> > <prefix>.paf")
        
#         sys.exit(1)
    
#     prefix = "prefix"
#     ref_file = bed_eval = raw_paf_file = corr_paf_file = output = None
#     for opt, arg in opts:
#         if opt in ['-o','--prefix']: prefix = arg
#         elif opt in ['-h','--hp']: ref_file = arg
#         elif opt in ['-r','--raw']: raw_paf_file = arg
#         elif opt in ['-c','--corrected']: corr_paf_file = arg
#         elif opt in ['-b','--specbed']: bed_eval = True    
    
# if __name__ == "__main__":
#     import sys
#     import getopt
#     import re
#     import matplotlib.pyplot as plt
#     from io import StringIO
#     from itertools import chain
#     from collections import defaultdict
#     import timeit
    
#     starttime = timeit.default_timer()
#     print("The start time is :",starttime)
    
#     main(sys.argv)
    
#     endtime = timeit.default_timer()
#     print("The end time is :",endtime)
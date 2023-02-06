#!/usr/bin/env python

# Global var that stores the most detailed information of hifieval
summary = "\t".join(["readName", "raw_mapped_chr", "raw_start", "raw_end", "corrected_mapped_chr", "corrected_start", "corrected_end", "num_oc", "num_uc", "num_cc\n"])

def add_summary(*args):
    global summary
    summary += "\t".join(args)
    summary += "\n"

# ------1. Class Objects for Error Correction Evaluation------ 
class ReadError(object):
    """This internal class object ReadError parses .paf files."""
    def __init__(self, readName, strand, chrName, start, end):
        self.readName = readName
        self.strand   = strand
        self.chrName  = chrName
        self.start    = int(start)
        self.end      = int(end)
        self.inPos    = []
        self.delPos   = []
        self.misPos   = []
    
    def __str__(self):
        return '\t'.join([self.readName,self.strand,self.chrName])
    
    def get_error(self, read_cs):
        """Get chromosome positions of reads mis-alignment"""
        tmp = read_cs
        if tmp.startswith('cs:Z:'):
            tmp = tmp[5:]
        else:
            eprint("Error in read {}'s cs: {}".format(self.readName, tmp))
            return self
        num = 0
        while len(tmp) != 0:
            # identical
            if tmp[0] == ':':
                identical = re.search(r'\d+', tmp).group()
                num += int(identical)
                tmp = re.sub(r'.', '', tmp, count = 1+len(identical))
            # mismatch
            elif tmp[0] == '*':
                num += 1
                self.misPos.append(num+self.start)
                tmp = re.sub(r'.', '', tmp, count = 3)
            # insertion
            elif tmp[0] == "+":
                insertion = ''
                for ch in tmp[1:]:
                    if ch.isalpha():
                        insertion += ch
                    else: break
                num += 1
                self.inPos.append(num+self.start)
                tmp = re.sub(r'.', '', tmp, count = 1+len(insertion))
            # deletion
            elif tmp[0] == "-":
                deletion = ''
                for ch in tmp[1:]:
                    if ch.isalpha():
                        deletion += ch
                    else: break
                num += 1
                self.delPos.append(num+self.start)
                tmp = re.sub(r'.', '', tmp, count = 1+len(deletion))
            else:
                eprint("Unexpected error tag {} in read {}'s cs: {}".format(tmp[0], self.readName, read_cs))
                break
        return self
    
    def get_reverse_error(self,read_cs):
        read_len = int(self.end)-int(self.start)+1
        tmp = read_cs
        if tmp.startswith('cs:Z:'):
            tmp = tmp[5:]
        else:
            eprint("Error in read {}'s cs: {}".format(self.readName, tmp))
            return self
        num = 0
        while len(tmp) != 0:
            # identical
            if tmp[0] == ':':
                identical = re.search(r'\d+', tmp).group()
                num += int(identical)
                tmp = re.sub(r'.', '', tmp, count = 1+len(identical))
            # mismatch
            elif tmp[0] == '*':
                num += 1
                self.misPos.append(read_len-num+1+self.start)
                tmp = re.sub(r'.', '', tmp, count = 3)
            # insertion
            elif tmp[0] == "+":
                insertion = ''
                for ch in tmp[1:]:
                    if ch.isalpha():
                        insertion += ch
                    else: 
                        break
                num += 1
                self.inPos.append(read_len-num+1+self.start)
                tmp = re.sub(r'.', '', tmp, count = 1+len(insertion))
            # deletion
            elif tmp[0] == "-":
                deletion = ''
                for ch in tmp[1:]:
                    if ch.isalpha():
                        deletion += ch
                    else: 
                        break
                num += 1
                self.delPos.append(read_len-num+1+self.start)
                tmp = re.sub(r'.', '', tmp, count = 1+len(deletion))
            else:
                eprint("Unexpected error tag {} in read {}'s cs: {}".format(tmp[0], self.readName, read_cs))
                break
        return self
    
    def all_error(self):
        e = self.inPos + self.delPos + self.misPos
        return e
    
class CorrectionStat(object):
    """This class object CorrectionStat compares ReadError from raw and corrected."""
    def __init__(self, chrName, readNum=0):
        self.chrName = chrName
        self.readNum = readNum
        self.ocPos   = [] # over-correction
        self.ucPos   = [] # under-correction
        self.ccPos   = [] # correct-correction
    
    def __str__(self):
        o = 'Chromosome {} has {} reads corrected by EC tool.'.format(self.chrName, str(self.readNum))
        return o
    
    def add_pos(self, ocPos, ucPos, ccPos):
        self.ocPos.append(ocPos)
        self.ucPos.append(ucPos)
        self.ccPos.append(ccPos)
        self.readNum += 1
        return self
    
    def read2chr(self):
        """Once user no longer needs read-lvl positions, collapse them into chr-lvl"""
        # This is non-reversible. 
        # Should think of a better way to store read- and chr- level stat
        self.ocPos = list(chain.from_iterable(self.ocPos))
        self.ucPos = list(chain.from_iterable(self.ucPos))
        self.ccPos = list(chain.from_iterable(self.ccPos))
        return self
    
    def get_FDR(self):
        # FDR = FP / (TP+FP)
        return len(self.ocPos)/(len(self.ccPos)+len(self.ocPos))
    
    def get_FNR(self):
        # FNR = FN / (TP+FN)
        return len(self.ucPos)/(len(self.ccPos)+len(self.ucPos))
    
    def get_TPR(self):
        # TPR = TP / (TP+FN)
        return len(self.ccPos)/(len(self.ccPos)+len(self.ucPos))


# ------2. Main functions for whole-genome Evaluation------
# ------2.1 Analyze the .paf files from minimap2.    ------
def paf_pairs(paf_file1, paf_file2):
    """
    Generator v3 for function error_correction_eval()
    Return a pair of list with selected columns from .paf files generated 
    by minimap2 mapping raw reads to ref genome and corrected reads
    to ref genome.
    """
    corr_pafs = readNames = None
    
    with open(paf_file2,'r') as paf2:
        # save the corrected reads paf for fast searching in raw reads pafs
        corr_pafs = [parse_read_paf(line) for line in paf2.readlines()]
        readNames = dict((paf[0], i) for i,paf in enumerate(corr_pafs))
        
    with open(paf_file1,'r') as paf1:
        prev_r = None
        
        for raw_line in paf1.readlines():
            raw_paf = parse_read_paf(raw_line, prev_r)
            if raw_paf is not None:    
                if raw_paf[0] in readNames:
                    # get corresponding corrected reads paf
                    corr_paf = corr_pafs[readNames[raw_paf[0]]]
                    prev_r = raw_paf
                else:
                    eprint('Read {} has been removed by the EC tool.'.format(raw_paf[0]))
                    prev_r = raw_paf
                    continue
                yield raw_paf, corr_paf

def parse_read_paf(line, prev_paf=None):
    """
    A Helper method to parse reads paf files for function error_correction_eval(),
    If prev_paf is always set to None, then secondary alignment is also considered.
    """
    line = line.rstrip()
    tmp = line.split('\t')
    
    # take the first appearing read in the paf file / primary alignment
    if prev_paf is not None:
        if tmp[0] == prev_paf[0] or len(tmp) < 24:
            return None
    
    indices = [0,4,5,7,8,10,-1]
    # readName, strand, chrName, start, end, nbase, cs
    readpaf = [tmp[index] for index in indices]
    
    return readpaf

def error_corr_helper(raw_plist, corr_plist):
    """
    A Helper method to calculate over-correction, under-correction, and correct-correction.
    Plist is a list of positions of all types of errors (in/del/mis) in one read
    """
    # The errors that are in raw but not corrected reads
    cc = set(raw_plist).difference(corr_plist)
    # Under-correction if corrected reads still has errors shared with raw reads
    uc = set(raw_plist).difference(cc)
    # Over-corrections are new in corrected reads but not in raw reads
    oc = set(corr_plist).difference(raw_plist)
    
    return list(oc), list(uc), list(cc)

def error_correction_eval(raw_paf_file, corr_paf_file):
    """
    Assess the results of read correction from error correction tool
    using the pairwise alignment info from raw and corrected reads to reference seq.
    oc: over-correction; uc: under-correction; cc: correct-correction
    """
    # The error correction stats objects stored in dict value
    # for each chromosome stored as dict key
    output = dict()
    ercor_tmp = None
    
    for raw_paf, corr_paf in paf_pairs(raw_paf_file, corr_paf_file):
        
        raw_error = ReadError(*raw_paf[0:5])
        if raw_error.strand == '+':
            raw_error = raw_error.get_error(raw_paf[6])
        else:
            raw_error = raw_error.get_reverse_error(raw_paf[6])

        corr_error = ReadError(*corr_paf[0:5])
        if corr_error.strand == '+':
            corr_error = corr_error.get_error(corr_paf[6])
        else:
            corr_error = corr_error.get_reverse_error(corr_paf[6])
        # check if reads are the same from raw and corrected
        if raw_error.readName != corr_error.readName:
            eprint('***Read name of two paf not matched***')
            eprint(raw_error)
            eprint(corr_error)
            return None
        oc, uc, cc = error_corr_helper(raw_error.all_error(), corr_error.all_error())
        
        add_summary(raw_error.readName, 
                    raw_error.chrName, str(raw_error.start), str(raw_error.end), 
                    corr_error.chrName, str(corr_error.start), str(corr_error.end), 
                    str(len(oc)), str(len(uc)), str(len(cc)))
        
        # assuming both raw and corrected reads are mapped to the same chromosome
        key = corr_error.chrName
        if key in output:
            output[key].add_pos(oc, uc, cc)
        else:
            output[key] = CorrectionStat(corr_error.chrName).add_pos(oc, uc, cc)

    return output

# ------2.2 Getting holistic stats for EC evaluation ------
def get_eval(ercor_output, chrlvl = True):
    """
    Evaluation metrics of one EC Tool.
    Returns a tab-delimited file with four columns: chrName, FDR, FNR, TPR.
    Used for generating bar plot.
    """
    all_chr_eval = defaultdict(float)
    txt = ''
    txt += '\t'.join(['chrName','FDR','FNR','TPR'])
    txt += '\n'
    for key in ercor_output.keys():
        tmp = ercor_output[key].read2chr()
        all_chr_eval['oc'] += len(tmp.ocPos)
        all_chr_eval['uc'] += len(tmp.ucPos)
        all_chr_eval['cc'] += len(tmp.ccPos)
        
        if chrlvl:
            fdr = tmp.get_FDR()
            fnr = tmp.get_FNR()
            tpr = tmp.get_TPR()
            txt += '\t'.join([key, str(fdr), str(fnr), str(tpr)])
            txt += '\n'
    
    fdr = all_chr_eval['oc']/(all_chr_eval['cc']+all_chr_eval['oc'])
    fnr = all_chr_eval['uc']/(all_chr_eval['cc']+all_chr_eval['uc'])
    tpr = all_chr_eval['cc']/(all_chr_eval['cc']+all_chr_eval['uc'])
    txt += '\t'.join(['all', str(fdr), str(fnr), str(tpr)])
    txt += '\n'
    return txt

def get_readlvl_eval(ercor_output, chrlvl = True):
    """
    Count how many corrected reads have 1 oc/uc, 2 oc/uc, ...
    Used for generating a histogram figure
    """
    oc_cntall = defaultdict(int)
    uc_cntall = defaultdict(int)
    txt = ''
    for key in ercor_output.keys():
        oc_cnt = defaultdict(int)
        uc_cnt = defaultdict(int)
        for read_oc in ercor_output[key].ocPos:
            if len(read_oc) > 0:
                oc_cnt[len(read_oc)] += 1
                oc_cntall[len(read_oc)] += 1
        for read_uc in ercor_output[key].ucPos:
            if len(read_uc) > 0:
                uc_cnt[len(read_uc)] += 1
                uc_cntall[len(read_uc)] += 1
        
        if chrlvl:
            txt += key+'_oc\t'+'\t'.join([str(oc_cnt[i]) 
                                          for i in sorted(oc_cnt.keys())])+'\n'
            txt += key+'_uc\t'+'\t'.join([str(uc_cnt[i]) 
                                          for i in sorted(uc_cnt.keys())])+'\n'
    
    txt += 'all_oc\t'+'\t'.join([str(oc_cntall[i]) 
                                 for i in sorted(oc_cntall.keys())])+'\n'
    txt += 'all_uc\t'+'\t'.join([str(uc_cntall[i]) 
                                 for i in sorted(uc_cntall.keys())])+'\n'
    
    if len(oc_cntall.keys()) > len(uc_cntall.keys()):
        cnt_nums = list(sorted(oc_cntall.keys()))
    else:
        cnt_nums = list(sorted(uc_cntall.keys()))
    txt = 'chr\t'+'\t'.join(cnt_nums) + '\n' + txt
    return txt


# ------3. Get EC performance in homopolymer regions ------
class UnitSTR(object):
    "Object class for Short Tandem Repeat info for one chromosome"
    def __init__(self, chrName):
        self.chrName = chrName
        self.strPos = []
        # ***Note: Not sure if we should care about the unit nucleotides 
        
    def add_strPos(self, start, end):
        self.strPos.append((start, end))
        return self
    
    def __str__(self):
        out = self.chrName + "\nNumber of Homopolymers:" + str(len(self.strPos))
        return out
    
    def to_bed(self):
        out = ''
        for pos in self.strPos:
            out += '\t'.join([self.chrName, str(pos[0]),str(pos[1])])
            out += '\n'
        return out
        
def find_homopolymers(fa_file, bed_file):
    chrName = unit = start = end = -1
    l = 3
    # the temporary UnitSTR object that will be stored in hp_dict
    tmp = None
    # a dict that store homopolymers for each chromosome
    hp_dict = {}
    with open(bed_file, 'w') as hp_bed:
        with open(fa_file,'r') as f:
            for line in f:
                if line[0] == '>':
                    if end - start > l:
                        tmp = tmp.add_strPos(start, end)
                        
                    if chrName != -1:
                        hp_dict[chrName] = tmp
                        hp_bed.write(tmp.to_bed())
                        
                    chrName = line.strip().split()[0][1:]
                    tmp = UnitSTR(chrName)
                    unit = -1
                    start = end = 0
                else:
                    for b in line.strip():
                        if b != unit:
                            if end - start > l:
                                tmp = tmp.add_strPos(start, end)
                            start = end
                            unit = b
                        end += 1
            if end - start > l:
                tmp = tmp.add_strPos(start, end)
            hp_dict[chrName] = tmp
            hp_bed.write(tmp.to_bed())
    return hp_dict


def hp_error_chr_eval(correction_dict, hp_dict):
    """hp_error_eval for each chromosome"""
    # check if the two dicts contains the same chromosomes
    if correction_dict.keys() != hp_dict.keys():
        print("Chromosomes don't match for Homopolymer Evaluation.")
        return
    
    output = dict()
    
    for key in correction_dict.keys():
        # for each chromosome
        corr_err = sorted(correction_dict[key].ocPos + correction_dict[key].ucPos)
        hp_pos = hp_dict[key].strPos
        
        hp_len_counter = defaultdict(int)
        hp_err = defaultdict(float) # error rate in each HP region
        i = j = 0
        
        while i < len(corr_err) and j < len(hp_pos):
            
            hp_len = hp_pos[j][1]-hp_pos[j][0]
            
            if corr_err[i] > hp_pos[j][1]:
                hp_len_counter[hp_len] += 1
                j += 1
            elif corr_err[i] < hp_pos[j][0]:
                i += 1
            else:
                hp_err[hp_len] += 1
                i += 1
                hp_len_counter[hp_len] += 1
                j += 1
        if len(hp_len_counter.keys()) == 0:
            continue
        hp_len_range = list(hp_len_counter.keys())
        max_len = max(hp_len_range)
        min_len = min(hp_len_range)
        
#         for i in range(max_len, min_len, -1):
#             hp_err[i-1] += hp_err[i]
        
        for i in range(min_len, max_len+1):
            if hp_len_counter[i]==0:
                continue
            hp_err[i] = hp_err[i] / (hp_len_counter[i])
        
        output[key] = hp_err
    return output

def hp_error_eval(correction_dict, hp_dict):
    
    # check if the two dicts contains the same chromosomes
    if correction_dict.keys() != hp_dict.keys():
        print("Chromosomes don't match for Homopolymer Evaluation.")
        print(correction_dict.keys(),hp_dict.keys())
        return
    
    hp_len_counter = defaultdict(int)
    hp_oc_err = defaultdict(int) # error rate in each HP region
    hp_uc_err = defaultdict(int)
    
    for key in correction_dict.keys():
        # for each chromosome
#         corr_err = sorted(correction_dict[key].ocPos + correction_dict[key].ucPos)
        oc = sorted(correction_dict[key].ocPos)
        uc = sorted(correction_dict[key].ucPos)
#         print(len(oc))
#         print(len(uc))
        hp_pos = hp_dict[key].strPos
        
        i = j = 0
        
        while i < len(oc) and j < len(hp_pos):
            
            hp_len = hp_pos[j][1]-hp_pos[j][0]
            
            if oc[i] > hp_pos[j][1]:
                hp_len_counter[hp_len] += 1
                j += 1
            elif oc[i] < hp_pos[j][0]:
                i += 1
            else:
                hp_oc_err[hp_len] += 1
                i += 1
                hp_len_counter[hp_len] += 1
                j += 1
                
        i = j = 0
        
        while i < len(uc) and j < len(hp_pos):
            
            hp_len = hp_pos[j][1]-hp_pos[j][0]
            
            if uc[i] > hp_pos[j][1]:
#                 hp_len_counter[hp_len] += 1
                j += 1
            elif uc[i] < hp_pos[j][0]:
                i += 1
            else:
                hp_uc_err[hp_len] += 1
                i += 1
#                 hp_len_counter[hp_len] += 1
                j += 1

    hp_len_range = list(hp_len_counter.keys())
    max_len = max(hp_len_range)
    min_len = min(hp_len_range)

    for i in range(min_len, max_len+1):
        if hp_len_counter[i]==0:
            continue
        hp_oc_err[i] = hp_oc_err[i] / hp_len_counter[i]
        hp_uc_err[i] = hp_uc_err[i] / hp_len_counter[i]
        
#         output.append(hp_err)

    return hp_oc_err,hp_uc_err,hp_len_counter

# ------ * Miscellaneous functions ------
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    
def dict2evalbed(correction_dict, prefix):
    uc_bedfile = open(prefix+'.uc.bed', 'w')
    lines = []
    for key in correction_dict.keys():
        tmp = correction_dict[key]
        for uc in tmp.ucPos:
            lines.append('\t'.join([key, str(uc), str(uc+1)])+'\n')
    uc_bedfile.writelines(lines)
    # Closing file
    uc_bedfile.close()
    
    oc_bedfile = open(prefix+'.oc.bed', 'w')
    lines = []
    for key in correction_dict.keys():
        tmp = correction_dict[key]
        for oc in tmp.ocPos:
            lines.append('\t'.join([key, str(oc), str(oc+1)])+'\n')
    oc_bedfile.writelines(lines)
    oc_bedfile.close()
    
    cc_bedfile = open(prefix+'.cc.bed', 'w')
    lines = []
    for key in correction_dict.keys():
        tmp = correction_dict[key]
        for cc in tmp.ccPos:
            lines.append('\t'.join([key, str(cc), str(cc+1)])+'\n')
    cc_bedfile.writelines(lines)
    cc_bedfile.close()
    
def specbed2dict(bed_file, chrNames):
    output = dict()
    for name in chrNames:
        output[name] = UnitSTR(name)
    with open(bed_file,'r') as f:
        while True:
            try:
                line = next(f)
                if line[0:3] != 'chr':
                    continue
                line = line.rstrip()
                line = line.split('\t')
                if line[0] == 'chrY':
                    continue
                output[line[0]].strPos.append((int(line[1]),int(line[2])))
            except StopIteration:
                break
    return output

def evalbed2dict(ocfile, ucfile, ccfile):
    output = dict()
    with open(ocfile,'r') as ocf:
        while True:
            try:
                tmp = next(ocf).rstrip().split('\t')
#                 print(tmp)
                key = tmp[0]
                if key in output:
                    output[key].ocPos.append(tmp[1])
                else:
                    output[key] = CorrectionStat(key)
                    output[key].ocPos.append(tmp[1])
            except StopIteration:
                break
    with open(ucfile,'r') as ucf:
        while True:
            try:
                tmp = next(ucf).rstrip().split('\t')
                key = tmp[0]
                output[key].ucPos.append(tmp[1])
            except StopIteration:
                break
                
    with open(ccfile,'r') as ccf:
        while True:
            try:
                tmp = next(ccf).rstrip().split('\t')
                key = tmp[0]
                output[key].ccPos.append(tmp[1])
            except StopIteration:
                break
    return output

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def main(argv):
    opts, args = getopt.getopt(argv[1:],"o:h:br:c:", 
                               ["prefix=","hp=","specbed","raw=","corrected="])
    
    if len(opts) < 2:
        print("Usage: hifieval.py  [options]")
        print("Options:")
        print("  -o STR      Output File Prefix")
        print("  -h STR      FASTA file with reference genome for evaluation in homopolymer region")
        print("  -b STR      BED file with specified regions for evaluation")
        print("  -r STR      PAF file aligned between raw reads and reference genome")
        print("  -c STR      PAF file aligned between corrected reads and reference genome")
        print("Minimap2 command for generating PAF file:")
        print("  minimap2 -t32 -cx map-hifi --secondary=no --paf-no-hit --cs <reference genome file> <reads file> > <prefix>.paf")
        
        sys.exit(1)
    
    prefix = "prefix"
    ref_file = bed_eval = raw_paf_file = corr_paf_file = output = None
    for opt, arg in opts:
        if opt in ['-o','--prefix']: prefix = arg
        elif opt in ['-h','--hp']: ref_file = arg
        elif opt in ['-r','--raw']: raw_paf_file = arg
        elif opt in ['-c','--corrected']: corr_paf_file = arg
        elif opt in ['-b','--specbed']: bed_eval = True
    
    if raw_paf_file is not None:
        output = error_correction_eval(raw_paf_file, corr_paf_file)
        # read correction status
        for i in sorted(output.keys(),key=natural_keys):
            eprint(output[i])
        # most detailed summary
        with open(prefix+'.summary.tsv', 'w') as f:
            f.write(summary)
        # generating histogram of readlvl evaluation
        with open(prefix+'.rdlvl.eval.tsv', 'w') as f:
            f.write(get_readlvl_eval(output))
        # overall metrics
        with open(prefix+'.metric.eval.tsv', 'w') as f:
            f.write(get_eval(output))
        
        dict2evalbed(output, prefix)

    if bed_eval:
        output = evalbed2dict(*args)
        if len(args) != 3:
            print("BED evaluation module failed.")
        else:
            with open(prefix+'.regionaleval.metric.tsv', 'w') as f:
                f.write(get_eval(output))        

    if ref_file is not None:
        hp_info = find_homopolymers(ref_file)
        hp_oc_err,hp_uc_err,hp_len_counter = hp_error_eval(output, hp_info)
        plt.scatter(hp_uc_err.keys(),
                    list(hp_uc_err.values()),alpha=.5, c='green',label='UC')
        plt.scatter(hp_oc_err.keys(),
                    list(hp_oc_err.values()),alpha=.5, c='darkorange',label='OC')
        plt.legend()
        plt.xlabel('Homopolymer Length')
        plt.ylabel('Error Rate')
        plt.yscale("log")
        plt.savefig(prefix+'HP.ErrorRate.png')

if __name__ == "__main__":
    import sys
    import getopt
    import re
    import matplotlib.pyplot as plt
    from io import StringIO
    from itertools import chain
    from collections import defaultdict
    main(sys.argv)
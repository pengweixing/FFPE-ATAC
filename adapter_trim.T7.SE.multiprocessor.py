#################################################
#  File Name:pyadapter_trim.xpw.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Sat 31 Oct 2020 05:15:48 PM UTC
#################################################

import sys
import os
import re
import argparse
import gzip
import bz2 
from Bio.Seq import Seq
from fuzzywuzzy import fuzz
import Levenshtein
from multiprocessing import Pool, Manager
import pandas as pd

def fargv():
    parser = argparse.ArgumentParser(usage="python pyadapter_trim.FFPE4.py -i reads_R1.fastq.gz")
    parser.add_argument('-i',"--R1",help="the input of fastq ", required=True)
    parser.add_argument('-m',"--mismatch",help="the number of mismatch",type=int,default=3)
    parser.add_argument('-p',"--processor",help="the number of processors",type=int,default=4)
    parser.add_argument('-len_mid',"--length_for_complete",help="the number of mismatch",type=int,default=15)
    parser.add_argument('-oh',"--expand_for_overhang",help="expand the flank of the reads",type=int,default=0)
    parser.add_argument('-len_min',"--minimal_adaptor",help="the minimal length of the adaptor",type=int,default=15)
    parser.add_argument('-len_out',"--length_for_output",help="the reads length for output",type=int,default=19)
    parser.add_argument('-thr',"--threhold_for_Levenshtein",help="the identity between adaptor and sequence",type=int,default=80)
    parser.add_argument('-cut_ada',"--cut_adaptor",help="the adaptor length for trimming",type=int,default=24)
    parser.add_argument('-a',"--adaptor",help="the adaptor's sequence",type=str,default='GGAGAAGATGTGTATAAGAGACAG')
    args = parser.parse_args()
    return args

def fuzz_align(reads,adaptor,mismatch):

    if not isinstance(reads,str):
        reads = reads.decode()
    for idx, base in enumerate(reads):  # loop through equal size windows
        dist = 0
        reads_subset = reads[idx:idx+len(adaptor)]
        if len(reads_subset)<len(adaptor):
            break
        if reads_subset == adaptor:
            return idx,dist
            break
        else:
            dist = Levenshtein.distance(reads_subset,adaptor)
            if dist <= mismatch:  # find first then break
                return idx,dist
                break
def trim(each_buffer,adaptor_list,adaptor_list_rc,args):
    
    mismatch = args.mismatch
    adaptor = args.adaptor
    length_for_complete = args.length_for_complete
    expand_for_overhang = args.expand_for_overhang
    minimal_adaptor = args.minimal_adaptor
    length_for_output = args.length_for_output
    threhold_for_Levenshtein = args.threhold_for_Levenshtein
    seq_header_o = []
    seq_o = []
    qual_header_o = []
    qual_o = []

    for line in each_buffer:
        seq_header = line[0]
        seq = next(each_buffer)[0]
        qual_header = next(each_buffer)[0]
        qual = next(each_buffer)[0]
        for each_subset in adaptor_list:
            if len(each_subset) >= length_for_complete:
                if fuzz.partial_ratio(each_subset,seq)>threhold_for_Levenshtein:
                    hold = fuzz_align(seq,each_subset,mismatch)
                    if hold:
                        idx,dist = hold
                        if dist<=mismatch:
                            seq = seq[idx+len(each_subset):]
                            qual = qual[idx+len(each_subset):]
                            break
            elif len(each_subset) < length_for_complete and len(each_subset)>minimal_adaptor:
                line_front = seq[0:len(each_subset)+expand_for_overhang]
                #dist = Levenshtein.distance(line_front, each_subset)
                dist = fuzz.partial_ratio(each_subset,line_front)
                dist = len(each_subset)*(1-dist)
                if dist<=mismatch:
                    seq = seq[len(each_subset)-1+expand_for_overhang:]
                    qual = qual[len(each_subset)-1+expand_for_overhang:]
                    break
        for each_subset in adaptor_list_rc:
            if len(each_subset) >= length_for_complete:
                hold = fuzz_align(seq,each_subset,mismatch)
                if fuzz.partial_ratio(each_subset,seq)>threhold_for_Levenshtein:
                    if hold:
                        idx,dist = hold
                        if dist <= mismatch:
                            seq = seq[0:idx]
                            qual = qual[0:idx]
                            break
            elif len(each_subset) < length_for_complete and len(each_subset)>minimal_adaptor:
                line_back = seq[-len(each_subset)-expand_for_overhang:]
                #dist = Levenshtein.distance(line_back, each_subset)
                dist = fuzz.partial_ratio(each_subset,line_front)
                dist = len(each_subset)*(1-dist)
                if dist<=mismatch:
                    seq = seq[0:-len(each_subset)-expand_for_overhang]
                    qual = qual[0:-len(each_subset)-expand_for_overhang]
                    break
        if  isinstance(seq_header,str):  
            seq_header = seq_header.encode() + b'\n'
            seq_header_o.append(seq_header)
            seq = seq.encode() + b'\n'
            seq_o.append(seq)
            qual_header = qual_header.encode()  + b'\n'
            qual_header_o.append(qual_header)
            qual = qual.encode() + b'\n'
            qual_o.append(qual)
        else:
            seq_header = seq_header + '\n'
            seq_header_o.append(seq_header)
            seq = seq + '\n'
            seq_o.append(seq)
            qual_header = qual_header  + '\n'
            qual_header_o.append(qual_header)
            qual = qual + '\n'
            qual_o.append(qual)
    all_trimmed_reads = []

    for index in range(len(seq_header_o)):
        if len(seq_o[index])>length_for_output:
            all_trimmed_reads.append([seq_header_o[index],seq_o[index],qual_header_o[index],qual_o[index]])
    return(all_trimmed_reads)
def gen_adaptor(adaptor):

    adaptor_subset = []
    for i in range(len(adaptor)):
        adaptor_subset.append(adaptor[i:i+len(adaptor)])
    return adaptor_subset
def gen_adaptor_rc(adaptor):

    adaptor_subset = []
    for i in range(len(adaptor)):
        adaptor_subset.append(adaptor[0:i+1])
    return adaptor_subset[::-1]

def main(kwargs):

    args = kwargs
    fa_input = args.R1
    processor = args.processor
    append = fa_input.split('.')[-1]
    mychunksize = 500000
    if append == "fastq":
        data = pd.read_csv(fa_input,sep="\t",chunksize=mychunksize,header=None) 
        R1_out = re.sub(".fastq", ".trim.fastq", fa_input)
        r1_write = open(R1_out, 'w')
    elif append == "fq":
        R1_out = re.sub(".fq", ".trim.fq", fa_input)
        r1_write = open(R1_out, 'w')
        data = pd.read_csv(fa_input,sep="\t",chunksize=mychunksize,header=None) 
    elif append == "gz":
        R1_out = re.sub(".gz", ".trim.gz", fa_input)
        r1_write = gzip.open(R1_out, 'wb')
        data = pd.read_csv(fa_input,compression='gzip',sep="\t",chunksize=mychunksize,header=None) 
    elif append == "bz2":
        R1_out = re.sub(".bz2", ".trim.bz2", fa_input)
        r1_write = gzip.open(R1_out, 'wb')
        data = pd.read_csv(fa_input,compression='bz2',sep="\t",chunksize=mychunksize,header=None) 
    else:
        sys.exit("ERROR! The input file2 must be a .fastq or .fastq.gz")

    cut_adaptor = args.cut_adaptor
    adaptor = args.adaptor
    adaptor_rc = str(Seq(adaptor).reverse_complement())
    adaptor_list = gen_adaptor(adaptor[-cut_adaptor:])
    adaptor_list_rc = gen_adaptor_rc(adaptor_rc[0:cut_adaptor])

    Processor = processor
    for each_chunk in data:
        each_chunk = each_chunk.values.tolist()
        pool = Pool(Processor)  
        reads = mychunksize/4       
        block = int(reads/Processor)
        results = []
        for num in range(Processor):
            start = num*block*4
            end = (num+1)*block*4
            each_buffer = iter(each_chunk[start:end])
            results.append(pool.apply_async(trim,(each_buffer,adaptor_list,adaptor_list_rc,args)))
        for each in results:
            for each2 in each.get():
                r1_write.write(each2[0])
                r1_write.write(each2[1])
                r1_write.write(each2[2])
                r1_write.write(each2[3])
        pool.close()
        pool.join()
    r1_write.close()
if __name__ == "__main__":
    kwargs = fargv()   
    main(kwargs)

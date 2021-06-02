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
from multiprocessing import Pool
import pandas as pd

def fargv():
    parser = argparse.ArgumentParser(usage="python adapter_trim.T7.PE.mutiprocessor.py -r1 reads_R1.fastq.gz -r2 reads_R2.fastq.gz")
    parser.add_argument('-r1',"--R1",help="the R1 of file ", required=True)
    parser.add_argument('-r2',"--R2",help="the R2 of file ", required=True)
    parser.add_argument('-p',"--processor",help="the number of processors default=4",type=int,default=4)
    parser.add_argument('-m',"--mismatch",help="the number of mismatch",type=int,default=3)
    parser.add_argument('-len_mid',"--length_for_complete",help="the cutoff for a complete adaptor sequence ",type=int,default=15)
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

def fuzz_align_compare(s_seq,l_seq,mismatch):
    for i, base in enumerate(l_seq):  # loop through equal size windows
        l_subset = l_seq[i:i+len(s_seq)]
        dist = Levenshtein.distance(l_subset, s_seq)
        if dist <= mismatch:  # find first then break
            return i, dist
            break

def trim_core(seq_in,qual_in,adaptor_list,adaptor_list_rc,args):

    mismatch = args.mismatch
    adaptor = args.adaptor
    length_for_complete = args.length_for_complete
    expand_for_overhang = args.expand_for_overhang
    minimal_adaptor = args.minimal_adaptor
    length_for_output = args.length_for_output
    threhold_for_Levenshtein = args.threhold_for_Levenshtein
    seq = seq_in
    qual = qual_in
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
    return seq,qual

def output(mybuffer,r1_write,r2_write):
    length = len(mybuffer.buffer_head1)
    length_for_output = mybuffer.length_for_output
    for i in range(length):
        seqhead1 = mybuffer.buffer_head1[i]
        seqhead2 = mybuffer.buffer_head2[i]
        seq1_o = mybuffer.buffer_seq1_o[i]
        seq2_o = mybuffer.buffer_seq2_o[i]
        qualhead1 = mybuffer.buffer_qualhead1[i]
        qualhead2 = mybuffer.buffer_qualhead2[i]          
        qual1_o = mybuffer.buffer_qual1_o[i]
        qual2_o = mybuffer.buffer_qual2_o[i]
        if  isinstance(seqhead1,str):  
            seqhead1 = seqhead1.encode() + b'\n'
            seq1_o = seq1_o.encode() + b'\n'
            qualhead1 = qualhead1.encode()  + b'\n'
            qual1_o = qual1_o.encode() + b'\n'
            seqhead2 = seqhead2.encode() + b'\n'
            seq2_o = seq2_o.encode() + b'\n'
            qualhead2 = qualhead2.encode()  + b'\n'
            qual2_o = qual2_o.encode() + b'\n'
        else:
            seqhead1 = seqhead1 + '\n'
            seq1_o = seq1_o + '\n'
            qualhead1 = qualhead1  + '\n'
            qual1_o = qual1_o + '\n'
            seqhead2 = seqhead2 + '\n'
            seq2_o = seq2_o + '\n'
            qualhead2 = qualhead2  + '\n'
            qual2_o = qual2_o + '\n'
        if len(seq1_o)>length_for_output and len(seq2_o)>length_for_output:
            r1_write.write(seqhead1)
            r1_write.write(seq1_o)
            r1_write.write(qualhead1)
            r1_write.write(qual1_o)
            r2_write.write(seqhead2)
            r2_write.write(seq2_o)
            r2_write.write(qualhead2)
            r2_write.write(qual2_o)
  # count += 1
def atac_trim(seq1_o,qual1_o,seq2_o,qual2_o):
    seq1 = seq1_o
    seq2 = seq2_o
    qual1 = qual1_o
    qual2 = qual2_o
    n = 18
    j = 0
    k = 0
    seq2 = seq2
    seq1 = seq1
    qual1 = qual1_o
    qual2 = qual2_o
    rc_seq2 = str(Seq(seq2[0:n]).reverse_complement())
    idx = seq1.rfind(rc_seq2) # look for perfect match
    mismatch = 1
    if idx > 0:
        j = j+1  # 0 mismatchs
    elif mismatch>0:
        hold = fuzz_align_compare(rc_seq2,seq1,mismatch)  # else allow for mismatch
        if hold:
            idx,mis=hold
            if mis == 1:
                k=k+1  # 1 mismatch

    # trim reads if idx exist
    if idx > 0:
        seq1 = seq1[0:idx+n-1]
        seq2 = seq2[0:idx+n-1]
        qual1 = qual1[0:idx+n-1]
        qual2 = qual2[0:idx+n-1]

    return seq1,qual1,seq2,qual2

class buffer:
    pass

mybuffer = buffer()
class buffer2:
    pass

def subread(args_list):
    mylen = len(args_list.R1_lines)
    count=1 
    mybuffer.buffer_seq1_o = []
    mybuffer.buffer_seq2_o = []
    mybuffer.buffer_qual1_o = []
    mybuffer.buffer_qual2_o = []
    mybuffer.buffer_head1 = []
    mybuffer.buffer_head2 = []
    mybuffer.buffer_qualhead1 = []
    mybuffer.buffer_qualhead2 = []
    for index in range(mylen):

        if count%4 == 1:
            seqhead1 = args_list.R1_lines[count-1]
            seqhead2 = args_list.R2_lines[count-1]
           # print(seqhead1)
        elif count%4 == 2:
            seq1 = args_list.R1_lines[count-1]
            seq2 = args_list.R2_lines[count-1]
          #  print(seq1)
        elif count%4 == 3:
            qualhead1 = args_list.R1_lines[count-1]
            qualhead2 = args_list.R2_lines[count-1]
        elif count%4 == 0:
            qual1 = args_list.R1_lines[count-1]
            qual2 = args_list.R2_lines[count-1]
            seq1_o,qual1_o = trim_core(seq1,qual1,args_list.adaptor_list,args_list.adaptor_list_rc,args_list.args)
            seq2_o,qual2_o = trim_core(seq2,qual2,args_list.adaptor_list,args_list.adaptor_list_rc,args_list.args)
            seq1_o2,qual1_o2,seq2_o2,qual2_o2 = atac_trim(seq1_o,qual1_o,seq2_o,qual2_o)

            mybuffer.buffer_seq1_o.append(seq1_o2)
            mybuffer.buffer_seq2_o.append(seq2_o2)
            mybuffer.buffer_qual1_o.append(qual1_o2)
            mybuffer.buffer_qual2_o.append(qual2_o2)
            mybuffer.buffer_head1.append(seqhead1)
            mybuffer.buffer_head2.append(seqhead2)
            mybuffer.buffer_qualhead1.append(qualhead1)
            mybuffer.buffer_qualhead2.append(qualhead2)
        count += 1
    mybuffer.length_for_output = args_list.args.length_for_output
    return mybuffer

def trim(R1_chunk,R2_chunk,adaptor_list,adaptor_list_rc,args):

    mismatch = args.mismatch
    adaptor = args.adaptor
    length_for_complete = args.length_for_complete
    expand_for_overhang = args.expand_for_overhang
    minimal_adaptor = args.minimal_adaptor
    length_for_output = args.length_for_output
    threhold_for_Levenshtein = args.threhold_for_Levenshtein
    para_list = []
    buffer2.adaptor_list = adaptor_list
    buffer2.adaptor_list_rc = adaptor_list_rc
    buffer2.args = args
    args_list = buffer2()
    args_list.R1_lines = []
    args_list.R2_lines = []
    rlt = []
    for r1_line,r2_line in zip(R1_chunk,R2_chunk):
        args_list.R1_lines.append(r1_line[0])
        args_list.R2_lines.append(r2_line[0])
  
    rlt.append(subread(args_list))
    
    return rlt

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
    R1_input = args.R1
    R2_input = args.R2
    append = R1_input.split('.')[-1]
    mychunksize = 100000
    processor = args.processor
    if append == "fastq":
      #  R1_rds = open(R1_input,'r')
     #   R2_rds = open(R2_input,'r')
        data_R1 = pd.read_csv(R1_input,sep="\t",chunksize=mychunksize,header=None) 
        data_R2 = pd.read_csv(R2_input,sep="\t",chunksize=mychunksize,header=None) 
        R1_out = re.sub(".fastq", ".trim.fastq", R1_input)
        R2_out = re.sub(".fastq", ".trim.fastq", R2_input)
        r1_write = open(R1_out, 'w')
        r2_write = open(R2_out, 'w')
    elif append == "fq":
      #  R1_rds = open(R1_input,'r')
      #  R2_rds = open(R2_input,'r')
        data_R1 = pd.read_csv(R1_input,sep="\t",chunksize=mychunksize,header=None) 
        data_R2 = pd.read_csv(R2_input,sep="\t",chunksize=mychunksize,header=None) 
        R1_out = re.sub(".fq", ".trim.fq", R1_input)
        R2_out = re.sub(".fq", ".trim.fq", R2_input)
        r1_write = open(R1_out, 'w')
        r2_write = open(R1_out, 'w')
    elif append == "gz":
        # R1_rds = gzip.open(R1_input,'r')
        # R2_rds = gzip.open(R2_input,'r')
        data_R1 = pd.read_csv(R1_input,compression='gzip',sep="\t",chunksize=mychunksize,header=None) 
        data_R2 = pd.read_csv(R2_input,compression='gzip',sep="\t",chunksize=mychunksize,header=None) 
        R1_out = re.sub(".gz", ".trim.gz", R1_input)
        R2_out = re.sub(".gz", ".trim.gz", R2_input)
        r1_write = gzip.open(R1_out, 'wb')
        r2_write = gzip.open(R2_out, 'wb')
    elif append == "bz2":
        # R1_rds = gzip.open(R1_input,'r')
        # R2_rds = gzip.open(R2_input,'r')
        data_R1 = pd.read_csv(R1_input,compression='bz2',sep="\t",chunksize=mychunksize,header=None) 
        data_R2 = pd.read_csv(R2_input,compression='bz2',sep="\t",chunksize=mychunksize,header=None) 
        R1_out = re.sub(".bz2", ".trim.bz2", R1_input)
        R2_out = re.sub(".bz2", ".trim.bz2", R2_input)
        r1_write = gzip.open(R1_out, 'wb')
        r2_write = gzip.open(R2_out, 'wb')
    else:
        sys.exit("ERROR! The input file2 must be a .fastq or .fastq.gz")
    cut_adaptor = args.cut_adaptor
    adaptor = args.adaptor
    adaptor_rc = str(Seq(adaptor).reverse_complement())
    adaptor_list = gen_adaptor(adaptor[-cut_adaptor:])
    adaptor_list_rc = gen_adaptor_rc(adaptor_rc[0:cut_adaptor])

    Processor = processor

    for data_R1_chunk,data_R2_chunk in zip(data_R1,data_R2):
        data_R1_chunk = data_R1_chunk.values.tolist()
        data_R2_chunk = data_R2_chunk.values.tolist()
        pool = Pool(Processor)
        reads = mychunksize/4
        block = int(reads/Processor)
        results = []
        for num in range(Processor):
            start = num*block*4
            end = (num+1)*block*4
            data_R1_chunk_each = data_R1_chunk[start:end]
            data_R2_chunk_each = data_R2_chunk[start:end]
            results.append(pool.apply_async(trim,(data_R1_chunk_each,data_R2_chunk_each,adaptor_list,adaptor_list_rc,args)))
        for rlts in results:
            for rlt in rlts.get():
                output(rlt,r1_write,r2_write) 
        pool.close()
        pool.join()

if __name__ == "__main__":
    kwargs = fargv()   
    main(kwargs)

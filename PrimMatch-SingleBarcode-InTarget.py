#!/usr/bin/python2.7

'''
    Date: 2017-11-04

    NOTE: we have modify the rules of sequence query.
    Assuming that every fastq sequence has a primer seq, we are sure
    that after query for "3+mismatch" bases, there is no possibility 
    to find the primer seq anymore!
'''

import sys
import os
import pdb
import pysam
import gzip
import os.path
from collections import defaultdict

# global variable
MAXPRIMLEN = 50
READLEN = 150
BLEN = 6 # bufer length for search

class FqIterator(object):
    def __init__(self, inFile):
        if inFile[-2:] == "gz":
           self.source = gzip.open(inFile, "rb")
        else:
           self.source = open(inFile, "rb")
        self.eof = False

    def next(self):
        group = []
        for i in xrange(4):
            try:
                line = self.source.next()
            except StopIteration:
                self.eof = True
                return
            if i == 0:
                group.append( line.split()[0] )
            else:
                group.append( line.rstrip() )
        return group

    def close(self):
        self.source.close()
        
   def GetPrim(primerfile):
    ''' primer = [('CCCTTCCACATA','CAAGTGAGGTCCTCAAAT', 99), (), ...]
        #Chr,Chain,Start,End,Seq,Chr,Primer_direction,Primer_start,Length,Primer_end,Distance,Distance1,PR_Start,PR_End,......
        chr7,+,55086916,55086943,TGACTCCGTCCAGTATTGATCGGGAGAG,,......
    '''
    primer = []

    with open(primerfile, "r") as fp:
        while (1):
            line = fp.readline()
            if not line.startswith('#'):
                break

        for line in fp.readlines():
            llist = line.split(",")
            primer.append( (llist[4], llist[0], int(llist[2]), int(llist[3])) )

    return primer


def SeqIndex(primseq, loc1, kmer):
    ''' primseq: primer sequence [str]
        loc1: primseq index in the primlist [int]
        //loc2: primseq index in the primer pair tuple [int]
        i: kmer seq index in the primer sequence
        index = {'ATGCG':[(loc1, i),(), ...], ...}
    '''
    index = {}
    cutlen = len(primseq) - kmer + 1

    for i in xrange(cutlen):
        key = primseq[i:i+kmer]
        if key not in index:
            index[key] = [(loc1, i)]
        else:
            index[key].append((loc1, i))

    return index
    
def Index(primlist, kmer=8):
    ''' index = {'ATGCG':[(0,0),(0,1)], 'TGCGA':[(),(),(),...], ...}
    '''
    index = {}
    pcount = len(primlist)

    for pindex in xrange(pcount):
        #for i in [0, 1]:
        tmpindex = SeqIndex(primlist[pindex][0], pindex, kmer)
        for key in tmpindex:
            if key not in index:
                index[key] = tmpindex[key]
            else:
                for ele in tmpindex[key]:
                    index[key].append(ele)
    return index


def RevComp(sequence):
    base = {'A':'T','T':'A','G':'C','C':'G', 'N':'N'}
    compseq = ''

    for b in sequence:
        compseq += base[b]

    return compseq[::-1]


def MisCheck(str1, str2, mismatch):
    misnum = 0
    len1 = len(str1)
    len2 = len(str2)

    minlen = len1 if len1 < len2 else len2
    for i in xrange(minlen):
        if str1[i] != str2[i]:
            misnum += 1

        if misnum > mismatch:
            return False

    return True 
    
 def SeqQuery(seq, index, primlist, mismatch, kmer):
    '''CORE FUNCTION
          seq: AGAAATTTGCGGAGTAAGTTGCGCTGGGGCTTTCGGCGGCGGCGATTTCGCC
       primer:      TTTGCGGAGTAAG
                    |           |
                  pstart       pend
       condition as above:  i = 5; h[2] = 0
    '''
    cutlen = 2 * BLEN

    for i in xrange(cutlen):
        try:
            hitlist = index[seq[i:i+kmer]]
        except KeyError:
            continue

        for h in hitlist:
            if i - h[1] >=0:
                string, primer = seq[i-h[1]:], primlist[h[0]][0]
                chrom, start, end = primlist[h[0]][1], primlist[h[0]][2], primlist[h[0]][3]
            elif (i - h[1] < 0) and (i - h[1] >= -mismatch):
                string, primer = seq, primlist[h[0]][0]
                chrom, start, end = primlist[h[0]][1], primlist[h[0]][2], primlist[h[0]][3]
            else:
                continue

            if MisCheck(string, primer, mismatch):
                pstart = 0 if i-h[1] < 0 else i-h[1]
                pend = pstart + len(primer) - 1
                primer_inseq = seq[pstart:pend+1]
                return (primer, primer_inseq, chrom, start, end)
    return ()
    
def GetM(cigar):
    length = 0
    for i in cigar:
        if i[0] == 0:
            length += i[1]

    return length
    
def main():
    args = sys.argv
    # args[1]: amplicon.txt
    # args[2]: fastq file

    if len(args) != 3:
        sys.stderr.write('Usage: python PrimMatch-1.2.1.py <amplicon.txt> <*.bam>\n')
        sys.exit(1)

    bamFile = pysam.AlignmentFile(args[2], "rb")
    primer = GetPrim(args[1])
    index = Index(primer, 8)
    filedone = False
    primer_count, total = {}, 0
    primer_read_num = defaultdict()         #primer_read_num[primer_in_amplicon] = [total_reads, in_target_reads, offtarget_reads]

    for read in bamFile.fetch(until_eof = True):
        try:
                #readCigar_M(read.cigar)
            if read.is_read1 and read.is_reverse:
                total += 1
                loc = SeqQuery(RevComp(read.seq)[:MAXPRIMLEN], index, primer, 3, 8)
                if loc:
                    if loc[0] not in primer_read_num:
                        primer_read_num[loc[0]] = [0,0,0]

                    primer_read_num[loc[0]][0] += 1
                    if (loc[2] == read.reference_name) and ((abs(loc[3]-(read.pos+len(read.seq)))<10) or ((abs(loc[3]-(read.pos+GetM(read.cigar)))<10))):
                        primer_read_num[loc[0]][1] += 1
                    else:
                        primer_read_num[loc[0]][2] += 1
            if read.is_read1 and (not read.is_reverse):
                total += 1
                loc = SeqQuery(read.seq[:MAXPRIMLEN], index, primer, 3, 8)
                if loc:
                    if loc[0] not in primer_read_num:
                        primer_read_num[loc[0]] = [0,0,0]
                    primer_read_num[loc[0]][0] += 1
                    if (loc[2] == read.reference_name) and (abs(loc[3]-read.pos)<10):
                        primer_read_num[loc[0]][1] += 1
                    else:
                        primer_read_num[loc[0]][2] += 1
                    '''
                    if loc[0] not in primer_count:
                        primer_count[loc[0]] = []
                        pos = loc[2] + "-" + str(loc[3]) + "-" + str(loc[4])    #pos in amplicon "chrom-start-end"
                        pos_in_read = read.reference_name + "-" + str(read.pos)
                        if [loc[1], pos_in_read] not in primer_count[loc[0]]:
                            primer_count[loc[0]].append([loc[1], pos_in_read])
                    else:
                        pos = loc[2] + "-" + str(loc[3]) + "-" + str(loc[4])
                        pos_in_read = read.reference_name + "-" + str(read.pos)
                        if [loc[1], pos_in_read] not in primer_count[loc[0]]:
                            primer_count[loc[0]].append([loc[1], pos_in_read])
                    '''
        except ValueError:
            continue

    '''output:primer_in_amplicon    primer_in_read  primer_in_read_pos
       TGAGTTCCTCAAAAGAGAAATCACGCATTTATGTT     TGAGTTCCTCAAAAGAGAAATCACGCATTTATGTT     chr8-98439990       
    '''
    print("primer\ttotal\tintarget\tofftarget")
    total_read1, total_intartget = 0, 0
    for key,value in primer_read_num.items():
        total_read1 += value[0]
        total_intartget +=  value[1]
        print("%s\t%d\t%d\t%d" %(key, value[0], value[1],value[2]))
                #print("%s\t%s\t%s\t%d\t%d" %(loc[0], loc[1], loc[2], loc[3], loc[4]))
    print("total reads in read1:%d" %(total))
    print("total find primer reads in read1:%d" %(total_read1))
    print("in target percent:%.2f%%" %(total_intartget/float(total)*100))

if __name__ == '__main__':
    main()    

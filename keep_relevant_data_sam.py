import csv
import sys
import math
import json
import multiprocessing as mp
import os
import numpy as np

samfile = sys.argv[1]
threads = int(sys.argv[2])
#skip_digits = int(sys.argv[3])
outfile = sys.argv[3]
datafolder = sys.argv[4]
folder = sys.argv[5]
subset1 = sys.argv[6]
subset2 = sys.argv[7]

numbers = ['0','1','2','3','4','5','6','7','8','9']
'''
def probabilityQ (X, b, p):
    pNew = getProbQuality(p)
    if X == b:
        qp = 1 - pNew
    else:
        qp = pNew/3
    return qp

def getProbQuality (q):
    qNew = ord(q)  - 33
    p = 10 ** (-np.float128(qNew)/10)
    return p
'''

def probability(a,b,qa,qb):
    qnewa = ord(qa) - 33
    qnewb = ord(qb) - 33
    pa = 10 ** (-np.float128(qnewa)/10)
    pb = 10 ** (-np.float128(qnewb)/10)
    if a == b:
        proba = (1-pa)*(1-pb) + pa*pb/3
    else:
        proba = (1-pa)*pb/3+pa*(1-pb)/3+2*pa*pb/9
    return proba

def getOverlapscore(tempSeq1,tempScore1,tempSeq2,tempScore2):
    prob = 0
    expect = 0
    nt = ["A","T","G","C"]
    L = min(len(tempSeq1),len(tempSeq2)) #POSITIONS
    for i in range(0, L): #LOOP THROUGH ALL POSITIONS OF THE MATCHES
        if (tempSeq1[i] not in nt) or (tempSeq2[i] not in nt):
            probabilityBase = 0.25
        else:
            probabilityBase = probability(tempSeq1[i],tempSeq2[i],tempScore1[i],tempScore2[i])
            '''
            probabilityBase = 0
            for n in nt:
                sc = (probabilityQ(n,tempSeq1[i],tempScore1[i]) * probabilityQ(n,tempSeq2[i],tempScore2[i])) #tempSeq1 => sequence of read1; tempScore1 => scores of read1...
                probabilityBase = probabilityBase + sc
            '''
        prob = prob + np.log(probabilityBase)
        expect = expect + probabilityBase
    expect = int(expect)
    overlapScore = round(np.exp(prob) ** (1/L),3)
    return (overlapScore,expect)
'''
def getExpNumMatches(tempSeq1,tempScore1,tempSeq2,tempScore2):
    expect = 0
    nt = ["A","T","G","C"]
    L = min(len(tempSeq1),len(tempSeq2)) #POSITIONS
    for i in range(0, L): #LOOP THROUGH ALL POSITIONS OF THE MATCHES
        if (tempSeq1[i] not in nt) or (tempSeq2[i] not in nt):
            probabilityBase = 0
        else:
            probabilityBase = 0
            for n in nt:
                sc = (probabilityQ(n,tempSeq1[i],tempScore1[i]) * probabilityQ(n,tempSeq2[i],tempScore2[i])) #tempSeq1 => sequence of read1; tempScore1 => scores of read1...
                probabilityBase = probabilityBase + sc
        expect = expect + probabilityBase
    expect = int(expect)
    return expect
'''
def rev(seq):
    revseq = seq[::-1]
    return(revseq)

def revcomp(seq):
    revc = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    revs = rev(seq)
    revcompl = []
    for n in revs:
        if n in revc.keys():
            revcompl.append(revc[n])
        else:
            revcompl.append(n)
    revcompl = ''.join(revcompl)
    return revcompl

def split_cigar(cig):
    cig_split = []
    current_split = []
    for c in range(len(cig)):
        current_split.append(cig[c])
        if cig[c] not in numbers:
            cig_split.append(''.join(current_split))
            current_split = []
    return cig_split


def test_extends_to_end(cig_split):
    count_begin = 0
    for i in range(len(cig_split)):
        if cig_split[i][-1] != 'M':
            count_begin += int(cig_split[i][:-1])
        else:
            break
    count_end = 0
    for i in range(len(cig_split)-1,-1,-1):
        if cig_split[i][-1] != 'M':
            count_end += int(cig_split[i][:-1])
        else:
            break
    if count_begin < 3 or count_end < 3:
        extends = 1
    else:
        extends = 0
    if count_begin <= count_end:
        start = 'begin'
    else:
        start = 'end'
    return (extends,start)

def find_ovlplen(cig_split,startdir):
    match_frac = []
    ovlplens = []
    ovlplen = 0
    match_bases = 0
    startpos = 0
    # Crop ends if they're not overlaps
    for i in range(len(cig_split)):
        if cig_split[i][-1] == 'M':
            start = i
            break
        else:
            if startdir == 'begin':
                startpos += int(cig_split[i][:-1])
    for i in range(len(cig_split)-1,-1,-1):
        if cig_split[i][-1] == 'M':
            end = i
            break
        else:
            if startdir == 'end':
                startpos += int(cig_split[i][:-1])
    cig_split_cr = cig_split[start:end+1]
    # Reverse CIGAR if the mapping is at the end of the query sequence
    if startdir == 'end':
        cig_split_cr = cig_split_cr[::-1]
    # Find the overlapping fragment
    for i in range(len(cig_split_cr)):
        ovlplen += int(cig_split_cr[i][:-1])
        if cig_split_cr[i][-1] == 'M':
            match_bases += int(cig_split_cr[i][:-1])
        ovlplens.append(ovlplen)
        match_frac.append(match_bases/ovlplen)
    ovlplen = 0
    for i in range(len(match_frac)-1,-1,-1):
        if match_frac[i] >= 0.80:
            ovlplen = ovlplens[i]
            break
    # Now find the start and end of the overlap on the query read
    endpos = startpos + ovlplen
    if startdir == 'end':
        readlen = sum(int(cig_split[i][:-1]) for i in range(len(cig_split)))
        endpos_temp = endpos
        endpos = readlen - startpos
        startpos = readlen - endpos_temp
    return (ovlplen,match_frac[i],startpos,endpos)

def processFile(fileID):
    with open('tempsam_'+subset1+'_'+subset2+'split'+str(fileID).zfill(2)+'.sam','r') as fr:
        with open('tempminimized_'+subset1+'_'+subset2+'_'+str(fileID).zfill(2)+'.sam','w') as fw:
            freader = csv.reader(fr, delimiter='\t')
            i = 0
            for row in freader:
                if row[2] != '*':
                    cig_split_ = split_cigar(row[5])
                    extends_, start_ = test_extends_to_end(cig_split_)
                else:
                    extends_ = 0
                if extends_ == 1:
                    query_seq = row[0]
                    ref_seq = row[2]
                    ovlplen_, match_, startpos_query, endpos_query = find_ovlplen(cig_split_,start_)
                    startpos_ref = int(row[3])-1
                    endpos_ref = startpos_ref + ovlplen_
                    if endpos_ref <= len(read_seq_phred[ref_seq][0]) and endpos_query <= len(read_seq_phred[query_seq][0]):
                        if bin(int(row[1]))[-5] == '0':
                            if (startpos_query < 3 and len(read_seq_phred[ref_seq][0]) - endpos_ref < 3) or (len(read_seq_phred[query_seq][0]) - endpos_query < 3 and startpos_ref < 3):
                                (overlapscore,expectation) = getOverlapscore(read_seq_phred[query_seq][0][startpos_query:endpos_query],read_seq_phred[query_seq][1][startpos_query:endpos_query],read_seq_phred[ref_seq][0][startpos_ref:endpos_ref],read_seq_phred[ref_seq][1][startpos_ref:endpos_ref])
 #                               expectation = getExpNumMatches(read_seq_phred[query_seq][0][startpos_query:endpos_query],read_seq_phred[query_seq][1][startpos_query:endpos_query],read_seq_phred[ref_seq][0][startpos_ref:endpos_ref],read_seq_phred[ref_seq][1][startpos_ref:endpos_ref])
                                fw.write(query_seq+'\t'+ref_seq+'\t'+str(ovlplen_)+'\t'+str(match_)+'\t'+str(round(overlapscore,3))+'\t'+str(round(expectation,3))+'\n') #+'\t'+read_seq_phred[query_seq][0][startpos_query:endpos_query]+'\t'+read_seq_phred[query_seq][1][startpos_query:endpos_query]+'\t'+read_seq_phred[ref_seq][0][startpos_ref:endpos_ref]+'\t'+read_seq_phred[ref_seq][1][startpos_ref:endpos_ref]+'\n')
                        else:
                            if (startpos_query < 3 and len(read_seq_phred[ref_seq][0]) - endpos_ref < 3) or (len(read_seq_phred[query_seq][0]) - endpos_query < 3 and startpos_ref < 3):
                                (overlapscore,expectation) = getOverlapscore(revcomp(read_seq_phred[query_seq][0])[startpos_query:endpos_query],rev(read_seq_phred[query_seq][1])[startpos_query:endpos_query],read_seq_phred[ref_seq][0][startpos_ref:endpos_ref],read_seq_phred[ref_seq][1][startpos_ref:endpos_ref])
#                                expectation = getExpNumMatches(revcomp(read_seq_phred[query_seq][0])[startpos_query:endpos_query],rev(read_seq_phred[query_seq][1])[startpos_query:endpos_query],read_seq_phred[ref_seq][0][startpos_ref:endpos_ref],read_seq_phred[ref_seq][1][startpos_ref:endpos_ref])
                                fw.write(query_seq+'\t'+ref_seq+'\t'+str(ovlplen_)+'\t'+str(match_)+'\t'+str(round(overlapscore,3))+'\t'+str(round(expectation,3))+'\n') #+'\t'+revcomp(read_seq_phred[query_seq][0])[startpos_query:endpos_query]+'\t'+rev(read_seq_phred[query_seq][1])[startpos_query:endpos_query]+'\t'+read_seq_phred[ref_seq][0][startpos_ref:endpos_ref]+'\t'+read_seq_phred[ref_seq][1][startpos_ref:endpos_ref]+'\n')
    return


if __name__ == '__main__':
    with open(samfile,'r') as f:
        num_lines = sum(1 for line in f)
    if num_lines > 0:
        read_seq_phred = json.load(open(folder+'read_seq_phred_'+subset1+'.json','r'))
        read_seq_phred2 = json.load(open(folder+'read_seq_phred_'+subset2+'.json','r'))
        read_seq_phred.update(read_seq_phred2)
        read_seq_phred2 = None
        num_lines_split = math.ceil(num_lines/threads)
        os.system('split -l '+str(num_lines_split)+' --numeric-suffixes=1 --additional-suffix=.sam '+folder+samfile+' tempsam_'+subset1+'_'+subset2+'split')
        po = mp.Pool(threads)
        for th in range(1,threads+1):
            po.apply_async(processFile, args=(th,))
        po.close()
        po.join()
        os.system('cat tempminimized_'+subset1+'_'+subset2+'* > '+folder+outfile)
        os.system('/bin/rm tempsam_'+subset1+'_'+subset2+'split*')
        os.system('/bin/rm tempminimized_'+subset1+'_'+subset2+'*')
    else:
        os.system('touch '+folder+outfile)

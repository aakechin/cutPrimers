# v13 - added checking of heterodimer formation of primers
# v14 - added ability to read and write to gzipped files; speed of processing was increased 10-times
# v15 - added ability to write untrimmed and trimmed reads to one file. Also added possibility that 3'-primer may be absent
# v16 - added ability to trim on the 3'-end only part of primer sequence
# v20 - added ability to trim primers in BAM-files
#       added ability to trim degenerate primers        

# Section of importing modules
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
bamNotAvailable=False
try:
    import pysam
except ModuleNotFoundError:
    print('WARNING: You do not have pysam module installed. You cannot trim primer sequences from BAM-files!')
    print('Install it or if you use Windows, you can not trim primers from BAM-files at all')
    bamNotAvailable=True
import glob,gzip,regex,time,argparse,math,hashlib,re
from multiprocessing import Pool,Queue
from multiprocessing.pool import ThreadPool
from itertools import repeat
from operator import itemgetter

def makeHashes(seq,k):
    # k is the length of parts
    subSeqs=[]
    h=[]
    lens=set([k])
    for i in range(len(seq)-k+1):
        h.append(hashlib.md5(seq[i:i+k].encode('utf-8')).hexdigest())
    return(h,lens)

def initializer(maxPrimerLen2,primerLocBuf2,errNumber2,primersR1_52,primersR1_32,primersR2_52,primersR2_32,
                primerR1_5_hashes2,primerR1_5_hashLens2,primerR2_5_hashes2,primerR2_5_hashLens2,
                primersFileR1_32,primersFileR2_52,primersFileR2_32,readsFileR22,primersStatistics2,
                idimer2,primer3absent2,minPrimer3Len2):
    global primersR1_5,primersR1_3,primersR2_5,primersR2_3,primersFileR1_3,primersFileR2_3,primersFileR2_5,readsFileR2
    global trimmedReadsR1,trimmedReadsR2,untrimmedReadsR1,untrimmedReadsR2
    global maxPrimerLen,q4,errNumber,primerLocBuf,readsPrimerNum,primersStatistics
    global primerR1_5_hashes,primerR2_5_hashes,primerR1_5_hashLens,primerR2_5_hashLens,primer3absent,idimer,minPrimer3Len
    maxPrimerLen=maxPrimerLen2
    primerLocBuf=primerLocBuf2
    errNumber=errNumber2
    primersR1_5=primersR1_52
    primersR1_3=primersR1_32
    primersR2_5=primersR2_52
    primersR2_3=primersR2_32
    primerR1_5_hashes=primerR1_5_hashes2; primerR1_5_hashLens=primerR1_5_hashLens2;
    primerR2_5_hashes=primerR2_5_hashes2; primerR2_5_hashLens=primerR2_5_hashLens2
    primersFileR1_3=primersFileR1_32
    primersFileR2_5=primersFileR2_52
    primersFileR2_3=primersFileR2_32
    readsFileR2=readsFileR22
    primersStatistics=primersStatistics2
    idimer=idimer2
    primer3absent=primer3absent2
    minPrimer3Len=minPrimer3Len2

# Section of functions
def showPercWork(done,allWork):
    percDoneWork=round((done/allWork)*100,2)
    sys.stdout.write("\r"+str(percDoneWork)+"%")
    sys.stdout.flush()

def revComplement(nuc):
    return(str(Seq(nuc).reverse_complement()))

def ambToRegList(seq):
    ambToNucs={'W':'[AT]',
               'Y':'[CT]',
               'R':'[AG]',
               'K':'[TG]',
               'S':'[CG]',
               'M':'[AC]',
               'B':'[TGC]',
               'H':'[ATC]',
               'D':'[AGT]',
               'V':'[AGC]',
               'N':'[ATGC]'}
    newSeq=[]
    for nuc in seq:
        if nuc in ambToNucs.keys():
            newSeq.append(ambToNucs[nuc])
        else:
            newSeq.append(nuc)
    return(''.join(newSeq))

def countDifs(s1,s2):
    a=pairwise2.align.globalms(s1,s2,2,-1,-1.53,0)
    maxSum=0
    k=0
    for i,b in enumerate(a):
        left=len(b[1])-len(b[1].lstrip('-'))+len(b[0])-len(b[0].lstrip('-'))
        right=len(b[1])-len(b[1].rstrip('-'))+len(b[0])-len(b[0].rstrip('-'))
        if left+right>maxSum:
            maxSum=left+right
            k=i
    ins=a[k][1].strip('-').count('-')
    dels=a[k][0].strip('-').count('-')
    left=max(len(a[k][1])-len(a[k][1].lstrip('-')),len(a[k][0])-len(a[k][0].lstrip('-')))
    right=max(len(a[k][1])-len(a[k][1].rstrip('-')),len(a[k][0])-len(a[k][0].rstrip('-')))
    if right==0:
        mism=sum(b!=c and c!='-' and b!='-' for b,c in zip(a[k][0][left:],a[k][1][left:]))
        return((mism,ins,dels,a[k][0][left:]))
    else:
        mism=sum(b!=c and c!='-' and b!='-' for b,c in zip(a[k][0][left:-right],a[k][1][left:-right]))
        return((mism,ins,dels,a[k][0][left:-right]))

def getErrors(s1,s2):
    # This function calculates number of errors between designed and sequenced primer sequences
    # s1 - initial sequence of primer
    # s2 - sequenced sequece of primer
    # Align them
    a=pairwise2.align.localms(s1,s2,2,-1,-1.53,0)
    maxSum=0
    k=0
    # First of all we detect the best alignment
    # and coordinates in range of which we will get mutations
    for i,b in enumerate(a):
        left=len(b[1])-len(b[1].lstrip('-'))+len(b[0])-len(b[0].lstrip('-'))
        right=len(b[1])-len(b[1].rstrip('-'))+len(b[0])-len(b[0].rstrip('-'))
        if left+right>maxSum:
            maxSum=left+right
            k=i
    poses=[] # poses - list of positions in sequences with mutations
    muts=[] # muts - mutations
    if right==0:
        s3=a[k][0][left:]
        s4=a[k][1][left:]
    else:
        s3=a[k][0][left:-right]
        s4=a[k][1][left:-right]
    for i,(b,c) in enumerate(zip(s3,s4)):
        if b!=c:
            poses.append(i+left+1)
            muts.append(b+'/'+c)
    return(poses,muts)

def trimPrimers(data):
    # This function get two records from both read files (R1 and R2)
    # and trim them
    # As a result it returns list
    #[trimmedReads,untrimmedReads]
    # resList is a variable with trimmed read sequences (0) and untrimmed read sequences (1)
    resList=[[None,None],[None,None]]
    r1,r2=data
    # Find primer at the 5'-end of R1 read
    readHashes=set()
    for l in primerR1_5_hashLens:
        hashes,lens=makeHashes(str(r1.seq[:maxPrimerLen+primerLocBuf]),l)
        readHashes.update(hashes)
    matchedPrimers={}
    for rh in readHashes:
        if rh in primerR1_5_hashes.keys():
            for a in primerR1_5_hashes[rh]:
                if a not in matchedPrimers.keys():
                    matchedPrimers[a]=1
                else:
                    matchedPrimers[a]+=1
    bestPrimer=None
    bestPrimerValue=None
    goodPrimers=[]
    goodPrimerNums=[]
    for key,item in sorted(matchedPrimers.items(),key=itemgetter(1),reverse=True):
        if bestPrimer==None:
            bestPrimer=key
            bestPrimerValue=item
            continue
        if item>=bestPrimerValue-1:
            goodPrimers.append(primersR1_5[key])
            goodPrimerNums.append(key)
        else: break
    if bestPrimer!=None:
        m1=regex.search(r''+primersR1_5[bestPrimer]+'{e<='+errNumber+'}',str(r1.seq[:maxPrimerLen+primerLocBuf]),flags=regex.BESTMATCH)
    else:
        return([[None,None],[r1,r2]],[],False)
##    m1=regex.search(r'(?:'+'|'.join(primersR1_5)+'){e<='+errNumber+'}',str(r1.seq[:maxPrimerLen+primerLocBuf]),flags=regex.BESTMATCH)
    # Use result of searching 5'-primer
    if m1==None:
        if len(goodPrimers)>0:
            m1=regex.search(r'(?:'+'|'.join(goodPrimers)+'){e<='+errNumber+'}',str(r1.seq[:maxPrimerLen+primerLocBuf]),flags=regex.BESTMATCH)
            if m1==None:
                # Save this pair of reads to untrimmed sequences
                return([[None,None],[r1,r2]],[],False)
            else:
                primerNum=goodPrimerNums[list(m1.groups()).index(m1[0])]
        else:
            return([[None,None],[r1,r2]],[],False)
    else:
        primerNum=bestPrimer
    # Find primer at the 5'-end of R2 read
    if primersFileR2_5:
        m3=regex.search(r'(?:'+primersR2_5[primerNum]+'){e<='+errNumber+'}',str(r2.seq[:maxPrimerLen+primerLocBuf]),flags=regex.BESTMATCH)
        if m3==None:
            # If user wants to identify hetero- and homodimers of primers
            if idimer:
                readHashes=set()
                for l in primerR2_5_hashLens:
                    hashes,lens=makeHashes(str(r2.seq[:maxPrimerLen+primerLocBuf]),l)
                    readHashes.update(hashes)
                matchedPrimers={}
                for rh in readHashes:
                    if rh in primerR2_5_hashes.keys():
                        for a in primerR2_5_hashes[rh]:
                            if a not in matchedPrimers.keys():
                                matchedPrimers[a]=1
                            else:
                                matchedPrimers[a]+=1
                bestPrimer=None
                bestPrimerValue=None
                goodPrimers=[]
                goodPrimerNums=[]
                for key,item in sorted(matchedPrimers.items(),key=itemgetter(1),reverse=True):
                    if bestPrimer==None:
                        bestPrimer=key
                        bestPrimerValue=item
                        continue
                    if item>=bestPrimerValue-1:
                        goodPrimers.append(primersR2_5[key])
                        goodPrimerNums.append(key)
                    else: break
                if bestPrimer!=None:
                    m3=regex.search(r''+primersR2_5[bestPrimer]+'{e<='+errNumber+'}',str(r2.seq[:maxPrimerLen+primerLocBuf]),flags=regex.BESTMATCH)
                else:
                    return([[None,None],[r1,r2]],[],False)
##                    m3=regex.search(r'(?:'+'|'.join(primersR2_5)+'){e<='+errNumber+'}',str(r2.seq[:maxPrimerLen+primerLocBuf]),flags=regex.BESTMATCH)
                # Use result of searching 5'-primer
                if m3==None:
                    if len(goodPrimers)>0:
                        m3=regex.search(r'(?:'+'|'.join(goodPrimers)+'){e<='+errNumber+'}',str(r2.seq[:maxPrimerLen+primerLocBuf]),flags=regex.BESTMATCH)
                        if m3==None:
                            # Save this pair of reads to untrimmed sequences
                            return([[None,None],[r1,r2]],[],False)
                        else:
                            primerNum2=goodPrimerNums[list(m3.groups()).index(m3[0])]
                    else:
                        return([[None,None],[r1,r2]],[],False)
                else:
                    primerNum2=bestPrimer
                # If we found two different 
                if primerNum!=primerNum2:
                    return([[None,None],[r1,r2]],[],[primerNum,primerNum2])
            else:
                # Save this pair of reads to untrimmed sequences
                return([[None,None],[r1,r2]],[],False)
        else:
            primerNum2=primerNum
    else:
        primerNum2=None
    # Find primer at the 3'-end of R1 read
    if primersFileR1_3:
        if not minPrimer3Len:
            m2=regex.search(r'(?:'+primersR1_3[primerNum]+'){e<='+errNumber+'}',str(r1.seq[-maxPrimerLen-primerLocBuf:]),flags=regex.BESTMATCH)
        else:
            errNumberDescreased=int(round(int(errNumber)*minPrimer3Len/len(primersR1_3[primerNum][:-2])))
            m2=regex.search(r'(?:'+primersR1_3[primerNum][:minPrimer3Len]+')){e<='+str(errNumberDescreased)+'}',str(r1.seq[-maxPrimerLen-primerLocBuf:]),flags=regex.BESTMATCH)
        if not primer3absent and m2==None:
            # Save this pair of reads to untrimmed sequences
            return([[None,None],[r1,r2]],[],[primerNum,primerNum2])
    # Find primer at the 3'-end of R2 read
    if primersFileR2_3:
        if not minPrimer3Len:
            m4=regex.search(r'(?:'+primersR2_3[primerNum]+'){e<='+errNumber+'}',str(r2.seq[-maxPrimerLen-primerLocBuf:]),flags=regex.BESTMATCH)
        else:
            errNumberDescreased=int(round(int(errNumber)*minPrimer3Len/len(primersR2_3[primerNum][:-2])))
            m4=regex.search(r'(?:'+primersR2_3[primerNum][:minPrimer3Len]+')){e<='+str(errNumberDescreased)+'}',str(r2.seq[-maxPrimerLen-primerLocBuf:]),flags=regex.BESTMATCH)
        if not primer3absent and m4==None:
            # Save this pair of reads to untrimmed sequences
            return([[None,None],[r1,r2]],[],[primerNum,primerNum2])
    # If all primers were found
    # Trim sequences of primers and write them to result file
    if primersFileR1_3 and m2!=None:
        r1.description+=':'+str(primerNum+1)
        resList[0][0]=r1[m1.span()[1]:len(r1.seq)-maxPrimerLen-primerLocBuf+m2.span()[0]]
    else:
        r1.description+=':'+str(primerNum+1)
        resList[0][0]=r1[m1.span()[1]:]
    if readsFileR2:
        if primersFileR2_3 and m4!=None:
            r2.description+=':'+str(primerNum+1)
            resList[0][1]=r2[m3.span()[1]:len(r2.seq)-maxPrimerLen-primerLocBuf+m4.span()[0]]
        elif primersFileR2_5:
            r2.description+=':'+str(primerNum+1)
            resList[0][1]=r2[m3.span()[1]:]
    # Save number of errors and primers sequences
    # [number of primer,difs1,difs2,difs3,difs4,]
    # Each dif is a set of (# of mismatches,# of insertions,# of deletions,primer_seq)
    if primersStatistics:
        difs1=countDifs(m1[0],primersR1_5[primerNum][1:-1])
        if primersFileR1_3 and m2!=None: difs2=countDifs(m2[0],primersR1_3[primerNum][1:-1])
        else: difs2=(0,0,0,'')
        if primersFileR2_5: difs3=countDifs(m3[0],primersR2_5[primerNum][1:-1])
        else: difs3=(0,0,0,'')
        if primersFileR2_3 and m4!=None: difs4=countDifs(m4[0],primersR2_3[primerNum][1:-1])
        else: difs4=(0,0,0,'')
        return (resList,[primerNum,difs1,difs2,difs3,difs4],False)
    else:
        return (resList,[],False)

def getPrimerNumFromRead(r,primersR1_5,primerR1_5_hashes,primerR1_5_hashLens,maxPrimerLen,primerLocBuf,readNum=1):
    readHashes=set()
    if readNum==1:
        for l in primerR1_5_hashLens:
            hashes,lens=makeHashes(str(r.seq[-maxPrimerLen-primerLocBuf:]),l)
            readHashes.update(hashes)
    else:
        for l in primerR1_5_hashLens:
            hashes,lens=makeHashes(str(r.seq[:maxPrimerLen+primerLocBuf]),l)
            readHashes.update(hashes)
    matchedPrimers={}
    for rh in readHashes:
        if rh in primerR1_5_hashes.keys():
            for a in primerR1_5_hashes[rh]:
                if a not in matchedPrimers.keys():
                    matchedPrimers[a]=1
                else:
                    matchedPrimers[a]+=1
    bestPrimer=None
    bestPrimerValue=None
    goodPrimers=[]
    goodPrimerNums=[]
    for key,item in sorted(matchedPrimers.items(),key=itemgetter(1),reverse=True):
        if bestPrimer==None:
            bestPrimer=key
            bestPrimerValue=item
            continue
        if item>=bestPrimerValue-1:
            goodPrimers.append(primersR1_5[key])
            goodPrimerNums.append(key)
        else: break
    if bestPrimer!=None:
        if readNum==1:
            m1=regex.search(r''+primersR1_5[bestPrimer]+'{e<='+errNumber+'}',str(r.seq[-maxPrimerLen-primerLocBuf:]),flags=regex.BESTMATCH)
        else:
            m1=regex.search(r''+primersR1_5[bestPrimer]+'{e<='+errNumber+'}',str(r.seq[:maxPrimerLen+primerLocBuf]),flags=regex.BESTMATCH)
    else:
        return(None,None)
    # Use result of searching 5'-primer
    if m1==None:
        if len(goodPrimers)>0:
            if readNum==1:
                m1=regex.search(r'(?:'+'|'.join(goodPrimers)+'){e<='+errNumber+'}',str(r.seq[-maxPrimerLen-primerLocBuf:]),flags=regex.BESTMATCH)
            else:
                m1=regex.search(r'(?:'+'|'.join(goodPrimers)+'){e<='+errNumber+'}',str(r.seq[:maxPrimerLen+primerLocBuf]),flags=regex.BESTMATCH)
            if m1==None:
                # Save this pair of reads to untrimmed sequences
                return(None,None)
            else:
                primerNum=goodPrimerNums[list(m1.groups()).index(m1[0])]
        else:
            return(None,None)
    else:
        primerNum=bestPrimer
    return(primerNum,m1)

def trimPrimersInBam(r,coordToPrimerNum,amplCoords,maxPrimerLen,primerLocBuf,errNumber,primersR1_5,primersR1_3,primersR2_5,primersR2_3,
                     primerR1_5_hashes,primerR1_5_hashLens,primerR2_5_hashes,primerR2_5_hashLens,
                     primersFileR1_3,primersFileR2_5,primersFileR2_3,primer3absent,minPrimer3Len,minReadLen,hardClipping):
    # This function get aligned read from BAM-file
    # As a result it returns read with newered information
    # Find primer at the 5'-end of R1 read
    if len(r.cigar)==0 or r.alen<minReadLen or 'H' in r.cigarstring:
        return(None,r)
    chrom=r.reference_name
    if chrom in coordToPrimerNum.keys():
        primerNumsCovered=getTheMostFitPrimerNumByPos(r.pos,r.alen,maxPrimerLen,coordToPrimerNum[chrom])
    amplBlockChrom=None
    amplBlockStart=None
    amplBlockEnd=None
    m1,m2,m3,m4=None,None,None,None
    # R1-read
    # If this read is R1 and it is reverse complement to reference
    # or other segment is not reverse complement to reference and it is R2
##    if ((str(bin(r.flag))[-7]=='1' and
##         str(bin(r.flag))[-5]=='1') or
##        (len(str(bin(r.flag)))>=10 and
##         str(bin(r.flag))[-8]=='1' and
##         str(bin(r.flag))[-6]=='0')):
    if str(bin(r.flag))[-7]=='1':
        if chrom not in coordToPrimerNum.keys() or len(primerNumsCovered)==0:
            primerNum,m1=getPrimerNumFromRead(r,primersR1_5,primerR1_5_hashes,primerR1_5_hashLens,maxPrimerLen,primerLocBuf,readNum=1)
            if primerNum is None:
                return(None,r)
        else:
            for primerNum,val in sorted(primerNumsCovered.items(),key=itemgetter(1),reverse=True):
                if str(bin(r.flag))[-5]=='1':
                    m1=regex.search(r''+primersR1_5[primerNum]+'{e<='+errNumber+'}',
                                    str(r.seq[-maxPrimerLen-primerLocBuf:]),flags=regex.BESTMATCH)
                else:
                    m1=regex.search(r''+primersR1_5[primerNum]+'{e<='+errNumber+'}',
                                    str(r.seq[:maxPrimerLen+primerLocBuf]),flags=regex.BESTMATCH)
                if m1!=None:
                    break
            if m1==None:
                return(None,r)
        amplBlockChrom,amplBlockStart,amplBlockEnd=amplCoords[primerNum]
        # Find primer at the 3'-end of R1 read
        if primersFileR1_3:
            if not minPrimer3Len:
                if str(bin(r.flag))[-5]=='1':
                    m2=regex.search(r'(?:'+primersR1_3[primerNum]+'){e<='+errNumber+'}',
                                    str(r.seq[:maxPrimerLen+primerLocBuf]),flags=regex.BESTMATCH)
                else:
                    m2=regex.search(r'(?:'+primersR1_3[primerNum]+'){e<='+errNumber+'}',
                                    str(r.seq[-maxPrimerLen-primerLocBuf:]),flags=regex.BESTMATCH)
            else:
                errNumberDescreased=int(round(int(errNumber)*minPrimer3Len/len(primersR1_3[primerNum][:-2])))
                if str(bin(r.flag))[-5]=='1':
                    m2=regex.search(r'(?:'+primersR1_3[primerNum][:minPrimer3Len]+')){e<='+str(errNumberDescreased)+'}',
                                    str(r.seq[:maxPrimerLen+primerLocBuf]),flags=regex.BESTMATCH)
                else:
                    m2=regex.search(r'(?:'+primersR1_3[primerNum][:minPrimer3Len]+')){e<='+str(errNumberDescreased)+'}',
                                    str(r.seq[-maxPrimerLen-primerLocBuf:]),flags=regex.BESTMATCH)
        # If all primers were found
        # Trim sequences of primers and write them to result file
        if m1!=None and primersFileR1_3 and m2!=None:
            # To perform soft clipping we should change cigar value and position
            if str(bin(r.flag))[-5]=='1':
                newRead=trimCigar(r,m2.span()[1],len(r.seq)-maxPrimerLen-primerLocBuf+m1.span()[0],
                                  amplBlockChrom,amplBlockStart,amplBlockEnd,hardClipping)
            else:
                newRead=trimCigar(r,m1.span()[1],len(r.seq)-maxPrimerLen-primerLocBuf+m2.span()[0],
                                  amplBlockChrom,amplBlockStart,amplBlockEnd,hardClipping)
        elif m1!=None and ((primersFileR1_3 and primer3absent) or
                           not primersFileR1_3):
            if str(bin(r.flag))[-5]=='1':
                newRead=trimCigar(r,None,len(r.seq)-maxPrimerLen-primerLocBuf+m1.span()[0],
                                  amplBlockChrom,amplBlockStart,amplBlockEnd,hardClipping)
            else:
                newRead=trimCigar(r,m1.span()[1],None,amplBlockChrom,amplBlockStart,amplBlockEnd,hardClipping)
        else:
            return(None,r)
    # R2-read
    # If other segment is reverse complement to reference and this read is R2
    # or it is R1 and it is not reverse complement to reference
##    elif ((len(str(bin(r.flag)))>=10 and
##           str(bin(r.flag))[-8]=='1' and
##           str(bin(r.flag))[-6]=='1') or
##          (str(bin(r.flag))[-7]=='1' and
##           str(bin(r.flag))[-5]=='0')):
    elif (len(str(bin(r.flag)))>=10 and
           str(bin(r.flag))[-8]=='1'):
        if chrom not in coordToPrimerNum.keys() or len(primerNumsCovered)==0:
            primerNum,m3=getPrimerNumFromRead(r,primersR2_5,primerR2_5_hashes,primerR2_5_hashLens,maxPrimerLen,primerLocBuf,readNum=2)
            if primerNum is None:
                return(None,r)
        else:
            for primerNum,val in sorted(primerNumsCovered.items(),key=itemgetter(1),reverse=True):
                if str(bin(r.flag))[-5]=='1':
                    m3=regex.search(r''+primersR2_5[primerNum]+'{e<='+errNumber+'}',
                                    str(r.seq[-maxPrimerLen-primerLocBuf:]),
                                    flags=regex.BESTMATCH)
                else:
                    m3=regex.search(r''+primersR2_5[primerNum]+'{e<='+errNumber+'}',
                                    str(r.seq[:maxPrimerLen+primerLocBuf]),
                                    flags=regex.BESTMATCH)
                if m3!=None:
                    break
            if m3==None:
                return(None,r)
        amplBlockChrom,amplBlockStart,amplBlockEnd=amplCoords[primerNum]
        # Find primer at the 3'-end of R2 read
        if primersFileR2_3:
            if not minPrimer3Len:
                if str(bin(r.flag))[-5]=='1':
                    m4=regex.search(r'(?:'+primersR2_3[primerNum]+'){e<='+errNumber+'}',
                                    str(r.seq[:maxPrimerLen+primerLocBuf]),
                                    flags=regex.BESTMATCH)
                else:
                    m4=regex.search(r'(?:'+primersR2_3[primerNum]+'){e<='+errNumber+'}',
                                    str(r.seq[-maxPrimerLen-primerLocBuf:]),
                                    flags=regex.BESTMATCH)
            else:
                errNumberDescreased=int(round(int(errNumber)*minPrimer3Len/len(primersR2_3[primerNum][:-2])))
                if str(bin(r.flag))[-5]=='1':
                    m4=regex.search(r'(?:'+primersR2_3[primerNum][:minPrimer3Len]+')){e<='+str(errNumberDescreased)+'}',
                                    str(r.seq[:maxPrimerLen+primerLocBuf]),
                                    flags=regex.BESTMATCH)
                else:
                    m4=regex.search(r'(?:'+primersR2_3[primerNum][:minPrimer3Len]+')){e<='+str(errNumberDescreased)+'}',
                                    str(r.seq[-maxPrimerLen-primerLocBuf:]),
                                    flags=regex.BESTMATCH)
        if m3!=None and primersFileR2_3 and m4!=None:
            if str(bin(r.flag))[-5]=='1':
                newRead=trimCigar(r,m4.span()[1],len(r.seq)-maxPrimerLen-primerLocBuf+m3.span()[0],
                                  amplBlockChrom,amplBlockStart,amplBlockEnd,hardClipping)
            else:
                newRead=trimCigar(r,m3.span()[1],len(r.seq)-maxPrimerLen-primerLocBuf+m4.span()[0],
                                  amplBlockChrom,amplBlockStart,amplBlockEnd,hardClipping)
        elif primersFileR2_5 and m3!=None and ((primersFileR2_3 and primer3absent) or
                                               not primersFileR2_3):
            if str(bin(r.flag))[-5]=='1':
                newRead=trimCigar(r,None,len(r.seq)-maxPrimerLen-primerLocBuf+m3.span()[0],
                                  amplBlockChrom,amplBlockStart,amplBlockEnd,hardClipping)
            else:
                newRead=trimCigar(r,m3.span()[1],None,
                                  amplBlockChrom,amplBlockStart,amplBlockEnd,hardClipping)
        elif primersFileR2_5:
            return(None,r)
    else:
        return(None,r)
    if newRead:
        if newRead.alen<minReadLen:
            return(None,r) 
        return(True,newRead)
    else:
        return(None,r)        

def getTheMostFitPrimerNumByPos(readStart,readLen,maxPrimerLen,coordToPrimerNumChrom):
##    poses=[]
    primerNumsCovered={}
    for pos in range(readStart+maxPrimerLen,readStart+readLen+1):
##    # Check the first base after maximal length of primer + position of read
##    poses.append(readStart+maxPrimerLen+1)
##    # Check the middle base between primer sequence and end of read
##    poses.append(int(round((readStart+maxPrimerLen+readStart+readLen)/2,0)))
##    # Check the middle base between two added positions
##    poses.append(int(round((poses[0]+poses[1])/2,0)))    
##    for pos in poses:
        if pos not in coordToPrimerNumChrom.keys():
            continue
        for primerNum in coordToPrimerNumChrom[pos]:
            if primerNum not in primerNumsCovered.keys():
                primerNumsCovered[primerNum]=1
            else:
                primerNumsCovered[primerNum]+=1
    return(primerNumsCovered)

def trimCigar(r,start=None,end=None,
              amplBlockChrom=None,
              amplBlockStart=None,amplBlockEnd=None,
              hardClipping=False):
    debug=None
##    debug='MN00909:41:000H2LNCW:1:21102:19422:4672'
    if ((start!=None and start<0) or
        (end!=None and end<0) or
        (start!=None and end!=None and start>=end)):
##        if debug!=None and r.qname==debug:
##            print('DEBUG INFO:')
##            print(r.qname)
##            print('For this read, start<0 or end<0 or start>=end')
##            print('Start:',start,'End:',end)
        return(None)
    oldCigar=r.cigar
    cigarStr=''
    amplLen=r.alen
    if hardClipping:
        clip='5'
        quals=r.qual
    else:
        clip='4'
    if start:
        newPos=r.pos+start
        amplLen-=start
    else:
        newPos=r.pos
    if not amplBlockStart or amplBlockChrom!=r.reference_name or r.pos+r.alen<amplBlockStart or r.pos>amplBlockEnd:
##        if debug!=None and r.qname==debug:
##            print('DEBUG INFO:')
##            print(r.qname)
##            print('For this read, amplBlockStart was not determined, amplBlockChrom!=reference name, '
##                  'r.pos+r.alen<amplBlockStart or r.pos>amplBlockEnd')
##            print('amplBlockStart:',amplBlockStart,'amplBlockChrom:',amplBlockChrom,'amplBlockEnd:',amplBlockEnd)
##            print('reference name:',r.reference_name,'r.pos:',r.pos,'r.alen:',r.alen)
        return(None)
        amplBlockStart=-1
    if not amplBlockEnd or amplBlockChrom!=r.reference_name or r.pos+r.alen<amplBlockStart or r.pos>amplBlockEnd:
##        if debug!=None and r.qname==debug:
##            print('DEBUG INFO:')
##            print(r.qname)
##            print('For this read, amplBlockEnd was not determined, amplBlockChrom!=reference name, '
##                  'r.pos+r.alen<amplBlockStart or r.pos>amplBlockEnd')
##            print('amplBlockStart:',amplBlockStart,'amplBlockChrom:',amplBlockChrom,'amplBlockEnd:',amplBlockEnd)
##            print('reference name:',r.reference_name,'r.pos:',r.pos,'r.alen:',r.alen)
        return(None)
        amplBlockEnd=10000000000
    for cig in oldCigar:
        cigarStr+=str(cig[0])*cig[1]
    # Remove hard clipped nucleotides
    cigarStr=cigarStr.strip('5')
    # Shift newPos onto number of S in cigar
    startSoftClipped=0
    try:
        if cigarStr[0]=='4' and (start or r.pos<amplBlockStart):
##            if debug!=None and r.qname==debug:
##                print('DEBUG INFO:')
##                print(r.qname)
##                print('For this read before shifting, cigarStr[0]==4 and (start or r.pos<amplBlockStart)')
##                print('cigarStr:',cigarStr,'start:',start,'r.pos:',r.pos,'amplBlockStart:',amplBlockStart)
##                print('startSoftClipped:',startSoftClipped,'newPos:',newPos,'amplLen:',amplLen)
            p=re.compile(cigarStr[0]+'+')
            m=p.findall(cigarStr)
            startSoftClipped=len(m[0])
            newPos-=len(m[0])
            amplLen+=len(m[0])
##            if debug!=None and r.qname==debug:
##                print('DEBUG INFO:')
##                print(r.qname)
##                print('For this read after shifting, cigarStr[0]==4 and (start or r.pos<amplBlockStart)')
##                print('cigarStr:',cigarStr,'start:',start,'r.pos:',r.pos,'amplBlockStart:',amplBlockStart)
##                print('startSoftClipped:',startSoftClipped,'newPos:',newPos,'amplLen:',amplLen)
    except IndexError:
        print('ERROR:',oldCigar)
        print(r)
        exit(0)
    # Shift newPos onto number of S in cigar
    endSoftClipped=0
    try:
        if cigarStr[-1]=='4' and (end or r.pos+r.alen>amplBlockEnd):
##            if debug!=None and r.qname==debug:
##                print('DEBUG INFO:')
##                print(r.qname)
##                print('For this read before shifting, cigarStr[-1]==4 and (end or r.pos+r.alen>amplBlockEnd)')
##                print('cigarStr:',cigarStr)
            p=re.compile(cigarStr[-1]+'+')
            m=p.findall(cigarStr)
            endSoftClipped=len(r.seq)-len(m[-1])
##            if debug!=None and r.qname==debug:
##                print('DEBUG INFO:')
##                print(r.qname)
##                print('For this read after shifting, cigarStr[-1]==4 and (end or r.pos+r.alen>amplBlockEnd)')
##                print('cigarStr:',cigarStr,'endSoftClipped',endSoftClipped)
    except IndexError:
        print('ERROR:',oldCigar)
        print(r)
        exit(0)
    # Check if soft-clipped sequence is alredy trimmed more than primer sequence
    if (amplBlockChrom==r.reference_name and
        r.pos>=amplBlockStart and
        (start==None or start<startSoftClipped)):
        start=None
        newPos=r.pos
##        if debug!=None and r.qname==debug:
##            print('DEBUG INFO:')
##            print(r.qname)
##            print('For this read, amplBlockChrom==r.reference_name and r.pos>=amplBlockStart and (start==None or start<startSoftClipped)')
##            print('start:',start,'newPos',newPos)
    if (amplBlockChrom==r.reference_name and
        r.pos+r.alen<=amplBlockEnd and
        (end==None or end<endSoftClipped)):
        end=None
##        if debug!=None and r.qname==debug:
##            print('DEBUG INFO:')
##            print(r.qname)
##            print('For this read, amplBlockChrom==r.reference_name and r.pos+r.alen<=amplBlockEnd and (end==None or end<endSoftClipped)')
##            print('end:',end)
    newCigarStr=''
    if (start or newPos<amplBlockStart-1) and (end or newPos+amplLen>amplBlockEnd):
##        if debug!=None and r.qname==debug:
##            print('DEBUG INFO:')
##            print(r.qname)
##            print('For this read, (start or newPos<amplBlockStart-1) and (end or newPos+amplLen>amplBlockEnd)')
##            print('start:',start,'end:',end,'newPos:',newPos,'amplLen:',amplLen,'amplBlockStart:',amplBlockStart,'amplBlockEnd:',amplBlockEnd)
        # We should consider deletions. We do not take them into account, while determining how we should trim cigar
        if not start:
            start=0
        # New start of matches in cigar 
        newStart=start
        if not end:
            end=len(cigarStr)
        newEnd=end
        startReady=False
        endReady=False
        if '2' in cigarStr or '1' in cigarStr:
##            if debug!=None and r.qname==debug:
##                print('DEBUG INFO:')
##                print(r.qname)
##                print('For this read, 2 in cigarStr or 1 in cigarStr')
##                print('cigarStr:',cigarStr)
            # Number of soft-clipped nucleotides before newStart
            ## It can differ due to deletions
            startSoftClippedNum=newStart
            endSoftClippedNum=newEnd
            for j in range(len(cigarStr)):
                if j>=newStart and newPos>=amplBlockStart-1:
                    startReady=True
                if not startReady:
                    if j<newStart:
                        # If there is deletion before end of primer
                        if cigarStr[j]=='2':
##                            if debug!=None and r.qname==debug:
##                                print('DEBUG INFO:')
##                                print(r.qname)
##                                print('For this read, when:')
##                                print('j:',j,'newStart:',newStart,'cigarStr[j]==2')
##                                print('newPos+=1','newStart+=1')
                            newPos+=1
                            newStart+=1
                        # If there is an insertion before end of primer
                        elif cigarStr[j]=='1':
##                            if debug!=None and r.qname==debug:
##                                print('DEBUG INFO:')
##                                print(r.qname)
##                                print('For this read, when:')
##                                print('j:',j,'newStart:',newStart,'cigarStr[j]==1')
##                                print('newPos-=1')
                            newPos-=1
                    # If primer sequence was removed, but read is wider than amplicon
                    ## -1 because when we read BAM-file with pysam, all positions starts from 0
                    elif newPos<amplBlockStart-1:
##                        if debug!=None and r.qname==debug:
##                            print('DEBUG INFO:')
##                            print(r.qname)
##                            print('For this read, when:')
##                            print('j:',j,'newStart:',newStart,'newPos<amplBlockStart-1')
##                            print('cigarStr[j]',cigarStr[j])
                        # Increase number of soft clipped nucleotides in all cases, except for deletions
                        if cigarStr[j]!='2':
                            startSoftClippedNum+=1
                        # Increase position in all cases, except for insertions
                        if cigarStr[j]!='1':
                            newPos+=1
                        newStart+=1
                # we do not use -1, because pos+length is the position that is next to the last position
                if j>newEnd and newPos-newStart+newEnd<=amplBlockEnd:
                    break
                if j<=newEnd:
                    # If there is deletion before start of primer
                    if cigarStr[j]=='2':# and newPos+newEnd<amplBlockEnd:
                        newEnd+=1
            if newPos-newStart+newEnd>amplBlockEnd:
                # startSoftClipped is a number of 'S' in the initial cigar
                newEnd=amplBlockEnd-r.pos+startSoftClipped
            endSoftClippedNum=len(cigarStr)-newEnd-cigarStr[newEnd:].count('2')
            newCigarStr=clip*startSoftClippedNum+cigarStr[newStart:newEnd]+clip*(endSoftClippedNum)
            if hardClipping:
                r.seq=r.seq[startSoftClippedNum:-endSoftClippedNum]
                r.qual=quals[startSoftClippedNum:-endSoftClippedNum]
        else:
            newStart=max(start,amplBlockStart-r.pos+startSoftClipped-1)
            newPos=max(newPos,amplBlockStart-1)
            newEnd=min(newEnd,amplBlockEnd-r.pos+startSoftClipped)
            if newStart>newEnd:
                newCigarStr=clip*len(cigarStr)
            else:
                newCigarStr=clip*newStart+cigarStr[newStart:newEnd]+clip*(len(cigarStr)-newEnd)
            if hardClipping:
                r.seq=r.seq[newStart:-(len(cigarStr)-newEnd)]
                r.qual=quals[newStart:-(len(cigarStr)-newEnd)]
            if startReady or endReady:
                print(startSoftClippedNum,endSoftClippedNum)
    elif (start or newPos<amplBlockStart-1):
##        if debug!=None and r.qname==debug:
##            print('DEBUG INFO:')
##            print(r.qname)
##            print('For this read, (start or newPos<amplBlockStart-1)')
##            print('start:',start,'end:',end,'newPos:',newPos,'amplLen:',amplLen,'amplBlockStart:',amplBlockStart,'amplBlockEnd:',amplBlockEnd)
        # We should consider deletions. We do not take them into account, while determining how we should trim cigar
        if not start:
            start=0
        # New start of matches in cigar 
        newStart=start
        # Number of soft-clipped nucleotides before newStart
        ## It can differ due to deletions
        startSoftClippedNum=newStart
        if '2' in cigarStr or '1' in cigarStr:
##            if debug!=None and r.qname==debug:
##                print('DEBUG INFO:')
##                print(r.qname)
##                print('For this read, 2 in cigarStr or 1 in cigarStr')
##                print('cigarStr:',cigarStr)
            for j in range(len(cigarStr)):
                if j>=newStart and newPos>=amplBlockStart-1:
                    break
                if j<newStart:
                    # If there is deletion before end of primer
                    if cigarStr[j]=='2':
##                        if debug!=None and r.qname==debug:
##                            print('DEBUG INFO:')
##                            print(r.qname)
##                            print('For this read, when:')
##                            print('j:',j,'newStart:',newStart,'newPos:',newPos)
##                            print('cigarStr[j]==2')
##                            print('newPos+=1','newStart+=1')
                        newPos+=1
                        newStart+=1
                    # If there is an insertion before end of primer
                    elif cigarStr[j]=='1':
##                        if debug!=None and r.qname==debug:
##                            print('DEBUG INFO:')
##                            print(r.qname)
##                            print('For this read, when:')
##                            print('j:',j,'newStart:',newStart,'newPos:',newPos)
##                            print('cigarStr[j]==1')
##                            print('newPos-=1')
                        newPos-=1
                # If primer sequence was removed, but read is wider than amplicon
                ## -1 because when we read BAM-file with pysam, all positions starts from 0
                elif newPos<amplBlockStart-1:
##                    if debug!=None and r.qname==debug:
##                        print('DEBUG INFO:')
##                        print(r.qname)
##                        print('For this read, when:')
##                        print('j:',j,'newStart:',newStart,'newPos:',newPos)
##                        print('newPos<amplBlockStart-1')
##                        print('cigarStr[j]',cigarStr[j])
                    # Increase number of soft clipped nucleotides in all cases, except for deletions
                    if cigarStr[j]!='2':
                        startSoftClippedNum+=1
                    # Increase position in all cases, except for insertions
                    if cigarStr[j]!='1':
                        newPos+=1
                    newStart+=1
##                    if debug!=None and r.qname==debug:
##                        print('DEBUG INFO:')
##                        print(r.qname)
##                        print('For this read, after all increases:')
##                        print('j:',j,'newStart:',newStart,'newPos:',newPos,'cigarStr[j]',cigarStr[j])
            newCigarStr=clip*startSoftClippedNum+cigarStr[newStart:]
            if hardClipping:
                r.seq=r.seq[startSoftClippedNum:]
                r.qual=quals[startSoftClippedNum:]
        else:
            newStart=max(start,amplBlockStart-r.pos-1)
            newPos=max(newPos,amplBlockStart-1)
            newCigarStr=clip*newStart+cigarStr[newStart:]
            if hardClipping:
                r.seq=r.seq[newStart:]
                r.qual=quals[newStart:]
    elif (end or newPos+amplLen>amplBlockEnd):
##        if debug!=None and r.qname==debug:
##            print('DEBUG INFO:')
##            print(r.qname)
##            print('For this read, (end or newPos+amplLen>amplBlockEnd)')
##            print('start:',start,'end:',end,'newPos:',newPos,'amplLen:',amplLen,'amplBlockStart:',amplBlockStart,'amplBlockEnd:',amplBlockEnd)
        # We should consider deletions. We do not count them when determine, how we should trim cigar
        delAfterEnd=0
        if not end:
            end=len(cigarStr)
        newEnd=end
        endSoftClippedNum=newEnd
        if '2' in cigarStr or '1' in cigarStr:
            for j in range(len(cigarStr)):
                # we do not use -1, because pos+length is the position that is next to the last position
                if j>newEnd and newPos+newEnd<=amplBlockEnd:
                    break
                if j<=newEnd:
                    # If there is deletion before start of primer
                    if cigarStr[j]=='2':# and newPos+newEnd<amplBlockEnd:
                        newEnd+=1
            if newPos+newEnd>amplBlockEnd:
                newEnd=amplBlockEnd-r.pos+startSoftClipped
            endSoftClippedNum=len(cigarStr)-newEnd-cigarStr[newEnd:].count('2')
            newCigarStr=cigarStr[:newEnd]+clip*(endSoftClippedNum)
            if hardClipping:
                r.seq=r.seq[:-endSoftClippedNum]
                r.qual=quals[:-endSoftClippedNum]
        else:
            newEnd=min(newEnd,amplBlockEnd-r.pos)
            newCigarStr=cigarStr[:newEnd]+clip*(len(cigarStr)-newEnd)
            if hardClipping:
                r.seq=r.seq[:-(len(cigarStr)-newEnd)]
                r.qual=quals[:-(len(cigarStr)-newEnd)]
    newCigar=[]
    if len(newCigarStr)==0 or set(newCigarStr)==set(['4']):
        return(None)
    while(len(newCigarStr)>0):
        p=re.compile(newCigarStr[0]+'+')
        m=p.findall(newCigarStr)
        newCigar.append((int(newCigarStr[0]),len(m[0])))
        newCigarStr=newCigarStr.replace(m[0],'',1)
    r.cigar=newCigar
    r.pos=newPos
    return(r)

def changeMateCoordinates(r,matePositions):
    if r.is_read1:
        if (r.query_name+'_2' in matePositions.keys() and
            matePositions[r.query_name+'_2'][0]!=None):
            r.mrnm=matePositions[r.query_name+'_2'][0]
            r.mpos=matePositions[r.query_name+'_2'][1]+1
        else:
            r.mate_is_unmapped=True
            r.mrnm=r.reference_id
            r.mpos=r.pos+1
    elif r.is_read2:
        if (r.query_name+'_1' in matePositions.keys() and
            matePositions[r.query_name+'_1'][0]!=None):
            r.mrnm=matePositions[r.query_name+'_1'][0]
            r.mpos=matePositions[r.query_name+'_1'][1]+1
        else:
            r.mate_is_unmapped=True
            r.mrnm=r.reference_id
            r.mpos=r.pos+1
    else:
        r.mate_is_unmapped=True
        r.mrnm=r.reference_id
        r.mpos=r.pos+1
    return(r)
        
    
if __name__ == "__main__":    
    # Section of reading arguments
    par=argparse.ArgumentParser(description='This script cuts primers from reads sequences')
    par.add_argument('--readsFile_r1','-r1',dest='readsFile1',type=str,help='file with R1 reads of one sample',required=False)
    par.add_argument('--readsFile_r2','-r2',dest='readsFile2',type=str,help='file with R2 reads of one sample',required=False)
    if not bamNotAvailable:
        par.add_argument('--bam-file','-bam',dest='bamFile',type=str,help='BAM-file in which you want to cut primers from reads',required=False)
        par.add_argument('--coordinates-file','-coord',dest='coordsFile',type=str,help='file with coordinates of amplicons in the BED-format (without column names and locations of primers): chromosome | start | end. '
                         'It is necessary for cutting primer sequences from BAM-file. Its order should be the same as for files with primer sequences',required=False)
        par.add_argument('--out-bam-file','-outbam',dest='outBamFile',type=str,help='name of file for output BAM-file with reads',required=False)
        par.add_argument('--out-untrimmed-bam-file','-outbam2',dest='outUntrimmedBamFile',type=str,help='name of file for output BAM-file with untrimmed reads. It is not required. If you do not use this parameter, all untrimmed reads will be lost',required=False)
        par.add_argument('--minimal-read-length','-minlen',dest='minReadLen',type=int,help='minimal length of read after trimming. Default: 10',default=10)
        par.add_argument('--hard-clipping','-hard',dest='hardClipping',action='store_true',help='use this parameter, if you want to trim reads wuith hard clipping. By default, primer sequences are trimmed with soft-clipping')
    par.add_argument('--primersFileR1_5','-pr15',dest='primersFileR1_5',type=str,help='fasta-file with sequences of primers on the 5\'-end of R1 reads',required=True)
    par.add_argument('--primersFileR2_5','-pr25',dest='primersFileR2_5',type=str,help='fasta-file with sequences of primers on the 5\'-end of R2 reads. Do not use this parameter if you have single-end reads',required=False)
    par.add_argument('--primersFileR1_3','-pr13',dest='primersFileR1_3',type=str,help='fasta-file with sequences of primers on the 3\'-end of R1 reads. It is not required. But if it is determined, -pr23 is necessary',required=False)
    par.add_argument('--primersFileR2_3','-pr23',dest='primersFileR2_3',type=str,help='fasta-file with sequences of primers on the 3\'-end of R2 reads',required=False)
    par.add_argument('--reads-with-forward','-freads',dest='forwardReadsNum',type=int,
                     help='number of file with reads that contain forward-primers (1 or 2). Default: 2',
                     default=2)
    par.add_argument('--trimmedReadsR1','-tr1',dest='trimmedReadsR1',type=str,help='name of file for trimmed R1 reads',required=False)
    par.add_argument('--trimmedReadsR2','-tr2',dest='trimmedReadsR2',type=str,help='name of file for trimmed R2 reads',required=False)
    par.add_argument('--untrimmedReadsR1','-utr1',dest='untrimmedReadsR1',type=str,help='name of file for untrimmed R1 reads. If you want to write reads that has not been trimmed to the same file as trimmed reads, type the same name',required=False)
    par.add_argument('--untrimmedReadsR2','-utr2',dest='untrimmedReadsR2',type=str,help='name of file for untrimmed R2 reads. If you want to write reads that has not been trimmed to the same file as trimmed reads, type the same name',required=False)
    par.add_argument('--primersStatistics','-stat',dest='primersStatistics',type=str,help='name of file for statistics of errors in primers. This works only for paired-end reads with primers at 3\'- and 5\'-ends',required=False)
    par.add_argument('--error-number','-err',dest='errNumber',type=int,help='number of errors (substitutions, insertions, deletions) that allowed during searching primer sequence in a read sequence. Default: 5',default=5)
    par.add_argument('--primer-location-buffer','-plb',dest='primerLocBuf',type=int,help='Buffer of primer location in the read from the start or end of read. If this value is zero, than cutPrimers will search for primer sequence in the region of the longest primer length. Default: 10',default=10)
    par.add_argument('--min-primer3-length','-primer3len',dest='minPrimer3Len',type=int,help="Minimal length of primer on the 3'-end to trim. Use this parameter, if you are ready to trim only part of primer sequence of the 3'-end of read")
    par.add_argument('--primer3-absent','-primer3',dest='primer3absent',action='store_true',help="if primer at the 3'-end may be absent, use this parameter")
    par.add_argument('--identify-dimers','-idimer',dest='idimer',type=str,help='use this parameter if you want to get statistics of homo- and heterodimer formation. Choose file to which statistics of primer-dimers will be written. This parameter may slightly decrease the speed of analysis')
    par.add_argument('--threads','-t',dest='threads',type=int,help='number of threads',default=2)
    args=par.parse_args()
    print('The command was:\n',' '.join(sys.argv))
    readsFileR1=args.readsFile1
    readsFileR2=args.readsFile2
    if not bamNotAvailable:
        bamFile=args.bamFile
        coordsFile=args.coordsFile
        minReadLen=args.minReadLen
    else:
        bamFile=None
    primersFileR1_5=args.primersFileR1_5
    primersFileR2_5=args.primersFileR2_5
    primersFileR1_3=args.primersFileR1_3
    primersFileR2_3=args.primersFileR2_3
    primer3absent=args.primer3absent
    minPrimer3Len=args.minPrimer3Len
    errNumber=str(args.errNumber)
    primerLocBuf=args.primerLocBuf
    primersStatistics=args.primersStatistics
    idimer=args.idimer
    if not bamFile and not readsFileR1:
        print('ERROR: you should use at least -r1 or -bam argument with native or aligned reads, respectively!')
        exit(1)
    if (primersFileR1_3 and not primersFileR1_5) or (not primersFileR2_5 and primersFileR2_3):
        print('ERROR: use of -pr13 or -pr23 should be accompanied by use of second one parameter for 5\'-end')
        exit(1)
    if not bamFile and ((not readsFileR2 and primersFileR2_5) or (not readsFileR2 and primersFileR2_3)):
        print('ERROR: use of -pr23 or -pr25 should be accompanied by use of readsFile2 parameter')
        exit(1)
    if readsFileR2 and not primersFileR2_5:
        print('ERROR: use of -r2 parameter should be accompanied by use of at least -pr25 parameter')
        exit(1)
    if readsFileR1 and (not args.trimmedReadsR1 or not args.untrimmedReadsR1):
        print('ERROR! If you use FASTQ-file as an input, you should use -tr1 and -utr1 arguments for output reads')
        exit(1)
    if bamFile and not args.outBamFile:
        print('ERROR! If you use BAM-file as an input, you should use -outbam for output aligned reads')
        exit(1)
    if bamFile and not coordsFile:
        print('ERROR! If you use BAM-file as an input, you should use -coord for file with coordinates of amplicons')
        exit(1)
    try:
        if args.trimmedReadsR1:
            if args.trimmedReadsR1[-3:]!='.gz':
                trimmedReadsR1=open(args.trimmedReadsR1,'w')
            else:
                trimmedReadsR1=gzip.open(args.trimmedReadsR1,'wt')
    except FileNotFoundError:
        print('########')
        print('ERROR! Could not create file:',args.trimmedReadsR1)
        print('########')
        exit(0)
    if args.untrimmedReadsR1 and args.trimmedReadsR1:
        if args.untrimmedReadsR1==args.trimmedReadsR1:
            untrimmedReadsR1=trimmedReadsR1
        else:
            try:
                if args.untrimmedReadsR1[-3:]!='.gz':
                    untrimmedReadsR1=open(args.untrimmedReadsR1,'w')
                else:
                    untrimmedReadsR1=gzip.open(args.untrimmedReadsR1,'wt')
            except FileNotFoundError:
                print('########')
                print('ERROR! Could not create file:',args.untrimmedReadsR1)
                print('########')
                exit(0)
    if args.trimmedReadsR2:
        try:
            if args.trimmedReadsR2[-3:]!='.gz':
                trimmedReadsR2=open(args.trimmedReadsR2,'w')
            else:
                trimmedReadsR2=gzip.open(args.trimmedReadsR2,'wt')
        except FileNotFoundError:
            print('########')
            print('ERROR! Could not create file:',args.trimmedReadsR2)
            print('########')
            exit(0)
    if args.untrimmedReadsR2:
        if args.untrimmedReadsR2==args.trimmedReadsR2:
            untrimmedReadsR2=trimmedReadsR2
        else:
            try:
                if args.untrimmedReadsR2[-3:]!='.gz':
                    untrimmedReadsR2=open(args.untrimmedReadsR2,'w')
                else:
                    untrimmedReadsR2=gzip.open(args.untrimmedReadsR2,'wt')
            except FileNotFoundError:
                print('########')
                print('ERROR! Could not create file:',args.untrimmedReadsR2)
                print('########')
                exit(0)
    if idimer and not readsFileR2:
        print('Warning! You did not provide R2-file so parameter "-idimer" will be ignored')
        idimer=None
    if idimer:
        try:
            idimerFile=open(idimer,'w')
        except FileNotFoundError:
            print('########')
            print('ERROR! Could not create file:',idimer)
            print('########')
            exit(0)
        primerDimers={}
    if primersStatistics:
        primersStatistics=open(args.primersStatistics,'w')
        primersStatisticsPos=open(args.primersStatistics[:-4]+'_poses.tab','w')
        primersStatisticsType=open(args.primersStatistics[:-4]+'_types.tab','w')
    threads=int(args.threads)

    # Read fasta-files with sequences of primers
    print('Reading files of primers...')
    lastPrimerNum=0
    # maxPrimerLen - variable that contains length of the longest primer
    maxPrimerLen=0
    # primers in R1 on the 5'-end
    primersR1_5=[]
    primersR1_5_names=[]
    primerR1_5_hashes={}
    primerR1_5_hashLens=set()
    primerR2_5_hashes={}
    primerR2_5_hashLens=set()
    i=0
    try:
        for r in SeqIO.parse(primersFileR1_5,'fasta'):
            primersR1_5_names.append(r.name)
            if bamFile and args.forwardReadsNum==2:
                primerSeq=str(r.seq.reverse_complement())
            else:
                primerSeq=str(r.seq)
            primersR1_5.append('('+ambToRegList(primerSeq)+')')
            partLens=math.floor(len(primerSeq)/(int(errNumber)+1))
            hashes,lens=makeHashes(primerSeq,partLens)
            primerR1_5_hashLens.update(lens)
            for h in hashes:
                if h in primerR1_5_hashes.keys():
                    primerR1_5_hashes[h].append(i)
                else:
                    primerR1_5_hashes[h]=[i]
            if len(r.seq)>maxPrimerLen:
                maxPrimerLen=len(r.seq)
            i+=1
    except FileNotFoundError:
        print('########')
        print('ERROR! File not found:',primersFileR1_5)
        print('########')
        exit(0)
    # primers in R2 on the 5'-end
    if primersFileR2_5:
        primersR2_5=[]
        primersR2_5_names=[]
        i=0
        try:
            for r in SeqIO.parse(primersFileR2_5,'fasta'):
                primersR2_5_names.append(r.name)
                if bamFile and args.forwardReadsNum==1:
                    primerSeq=str(r.seq.reverse_complement())
                else:
                    primerSeq=str(r.seq)
                primersR2_5.append('('+ambToRegList(primerSeq)+')')
                partLens=math.floor(len(primerSeq)/(int(errNumber)+1))
                hashes,lens=makeHashes(primerSeq,partLens)
                primerR2_5_hashLens.update(lens)
                for h in hashes:
                    if h in primerR2_5_hashes.keys():
                        primerR2_5_hashes[h].append(i)
                    else:
                        primerR2_5_hashes[h]=[i]
                if len(r.seq)>maxPrimerLen:
                    maxPrimerLen=len(r.seq)
                i+=1
        except FileNotFoundError:
            print('########')
            print('ERROR! File not found:',primersFileR2_5)
            print('########')
            exit(0)
    else:
        primersR2_5=None
    # primers in R1 on the 3'-end
    if primersFileR1_3:
        primersR1_3=[]
        primersR1_3_names=[]
        try:
            for r in SeqIO.parse(primersFileR1_3,'fasta'):
                primersR1_3_names.append(r.name)
                if bamFile and args.forwardReadsNum==2:
                    primerSeq=str(r.seq.reverse_complement())
                else:
                    primerSeq=str(r.seq)
                primersR1_3.append('('+ambToRegList(primerSeq)+')')
                if len(r.seq)>maxPrimerLen:
                    maxPrimerLen=len(r.seq)
        except FileNotFoundError:
            print('########')
            print('ERROR! File not found:',primersFileR1_3)
            print('########')
            exit(0)
    else:
        primersR1_3=None
    # primers in R2 on the 3'-end
    if primersFileR2_3:
        primersR2_3=[]
        primersR2_3_names=[]
        try:
            for r in SeqIO.parse(primersFileR2_3,'fasta'):
                primersR2_3_names.append(r.name)
                if bamFile and args.forwardReadsNum==1:
                    primerSeq=str(r.seq.reverse_complement())
                else:
                    primerSeq=str(r.seq)
                primersR2_3.append('('+ambToRegList(primerSeq)+')')
                if len(r.seq)>maxPrimerLen:
                    maxPrimerLen=len(r.seq)
        except FileNotFoundError:
            print('########')
            print('ERROR! File not found:',primersFileR2_3)
            print('########')
            exit(0)
    else:
        primersR2_3=None
    # Read file with R1 and R2 reads
    if readsFileR1:
        try:
            if readsFileR1[-3:]!='.gz':
                allWork=open(readsFileR1).read().count('\n')/4
            else:
                allWork=gzip.open(readsFileR1,'rt').read().count('\n')/4
        except FileNotFoundError:
            print('########')
            print('ERROR! Could not open file:',readsFileR1)
            print('########')
            exit(0)
    if readsFileR1:
        print('Reading input FASTQ-file(s)...')
        if readsFileR1[-3:]!='.gz':
            data1=SeqIO.parse(readsFileR1,'fastq')
        else:
            data1=SeqIO.parse(gzip.open(readsFileR1,'rt'),'fastq')
        if readsFileR2:
            try:
                if readsFileR2[-3:]!='.gz':
                    data2=SeqIO.parse(readsFileR2,'fastq')
                else:
                    data2=SeqIO.parse(gzip.open(readsFileR2,'rt'),'fastq')
            except FileNotFoundError:
                print('########')
                print('ERROR! Could not open file:',readsFileR2)
                print('########')
                exit(0)
        else:
            data2=['']*int(allWork)
    
        # Create Queue for storing result and Pool for multiprocessing
        primerErrorQ=[] 
        p=Pool(threads,initializer,(maxPrimerLen,primerLocBuf,errNumber,primersR1_5,primersR1_3,primersR2_5,primersR2_3,
                                    primerR1_5_hashes,primerR1_5_hashLens,primerR2_5_hashes,primerR2_5_hashLens,
                                    primersFileR1_3,primersFileR2_5,primersFileR2_3,readsFileR2,primersStatistics,
                                    idimer,primer3absent,minPrimer3Len))
        # Cutting primers and writing result immediately
        print('Trimming primers from reads...')
        doneWork=0
        showPercWork(0,allWork)
        for res in p.imap_unordered(trimPrimers,zip(data1,data2),10):
            doneWork+=1
            showPercWork(doneWork,allWork)
            if res[1]!=[]:
                primerErrorQ.append(res[1])
            if readsFileR2:
                if res[0][0][0] is not None and res[0][0][1] is not None:
                    SeqIO.write(res[0][0][0],trimmedReadsR1,'fastq')
                    SeqIO.write(res[0][0][1],trimmedReadsR2,'fastq')
                elif res[0][1][0] is not None and res[0][1][1] is not None:
                    # If user want to identify primer-dimers
                    if idimer and res[2]:
                        r1partSeq=str(res[0][1][0].seq[:40])
                        r2partSeq=revComplement(str(res[0][1][1].seq[:40]))
                        difs=countDifs(r1partSeq,r2partSeq)
                        if sum(difs[0:2])<=int(errNumber):
                            # and len(difs[3])>=len(primersR1_5[res[2][0]])
                            if primersR1_5_names[res[2][0]]+' & '+primersR2_5_names[res[2][1]] not in primerDimers.keys():
                                primerDimers[primersR1_5_names[res[2][0]]+' & '+primersR2_5_names[res[2][1]]]=1
                            else:
                                primerDimers[primersR1_5_names[res[2][0]]+' & '+primersR2_5_names[res[2][1]]]+=1
                        else:
                            SeqIO.write(res[0][1][0],untrimmedReadsR1,'fastq')
                            SeqIO.write(res[0][1][1],untrimmedReadsR2,'fastq')
                    else:
                        SeqIO.write(res[0][1][0],untrimmedReadsR1,'fastq')
                        SeqIO.write(res[0][1][1],untrimmedReadsR2,'fastq')
                            
                else:
                    print('ERROR: nor the 1st item of function result list or 2nd contains anything')
                    print(res)
                    exit(0)
            else:
                if res[0][0][0] is not None:
                    SeqIO.write(res[0][0][0],trimmedReadsR1,'fastq')
                elif res[0][1][0] is not None:
                    SeqIO.write(res[0][1][0],untrimmedReadsR1,'fastq')
                else:
                    print('ERROR: item of function result list contains anything')
                    print(res)
                    exit(0)
        print()
        # primersErrors is a dictionary that contains errors in primers
        if args.primersStatistics:
            primersErrors={}
            # primersErrorsPos is a dictionary that contains statistics about location
            # of errors
            primersErrorsPos=[{},{}]
            # primersErrorsType is a dictionary that contains statistics about type of error
            primersErrorsType=[{},{}]
            print('Counting errors...')
            for item in primerErrorQ:
                # If key for this primer has not been created, yet
                if not item[0] in primersErrors.keys():
                    # For each primer of each pair we will gather the following values:
                    # [(0)number of read pairs,
                    # (1)number of primers without errors,
                    # (2)number of primers with sequencing errors,
                    # (3)number of primers with synthesis errors
                    # The first item of list - F
                    # The second - R
                    primersErrors[item[0]]=[[0,0,0,0],[0,0,0,0]]
                    
    ##          R                           F_reverse_complement
    ## R1 5'---------________________________---------3'
    ## R2 5'---------________________________---------3'
    ##          F                           R_reverse_complement
                    
                # Increase number of read pairs
                primersErrors[item[0]][0][0]+=1
                primersErrors[item[0]][1][0]+=1
                # F-primers of pairs
                # The last variant is a case when we have single-end reads and 3' does not contain primer sequence
                if ((not primersFileR1_3 and primersFileR2_5 and item[3][0:3]==(0,0,0)) or
                    (primersFileR1_3 and primersFileR2_5 and item[3][0:3]==(0,0,0) and item[2][0:3]==(0,0,0)) or
                    (primersFileR1_5 and not primersFileR2_5 and not primersFileR1_3 and not primersFileR2_3)):
                    primersErrors[item[0]][0][1]+=1
                # If it was overlapping paired-end reads, we try to check if this is sequencing error
                elif primersFileR1_3 and primersFileR2_5 and primersFileR2_3 and item[2][3]!='' and item[3][3]!='':
                    # Reverse complement one of primer sequences
                    rev=str(Seq(item[2][3]).reverse_complement())
                    a=pairwise2.align.globalms(rev,item[3][3],2,-1,-1.53,-0.1)
                    # If found sequences are identical, it's a synthesis error
                    if list(a[0][0])==list(a[0][1]):
                        primersErrors[item[0]][0][3]+=1
                        # Now we want to save information about error's location
                        poses,muts=getErrors(primersR2_5[item[0]][1:-1],item[3][3])
                        for p in poses:
                            if p not in primersErrorsPos[1].keys():
                                primersErrorsPos[1][p]=1
                            else:
                                primersErrorsPos[1][p]+=1
                        for m in muts:
                            if m not in primersErrorsType[1].keys():
                                primersErrorsType[1][m]=1
                            else:
                                primersErrorsType[1][m]+=1
                    # Else it's a sequencing error
                    else:
                        primersErrors[item[0]][0][2]+=1
                # Else we just save it as sequencing error
                else:
                    primersErrors[item[0]][0][2]+=1
                # R-primers of pairs
                # For R-primer we always have sequence at least at 5' end of R1
                if ((not primersFileR2_3 and item[1][0:3]==(0,0,0)) or
                    (primersFileR2_3 and item[1][0:3]==(0,0,0) and item[4][0:3]==(0,0,0))):
                    primersErrors[item[0]][1][1]+=1
                # If it was overlapping paired-end reads, we try to check if this is sequencing error
                elif primersFileR1_3 and primersFileR2_5 and primersFileR2_3 and item[4][3]!='' and item[1][3]!='':
                    # Reverse complement one of primer sequences
                    rev=str(Seq(item[4][3]).reverse_complement())
                    a=pairwise2.align.globalms(rev,item[1][3],2,-1,-1.53,-0.1)
                    # If found sequences are identical, it's a synthesis error
                    try:
                        if list(a[0][0])==list(a[0][1]):
                            primersErrors[item[0]][1][3]+=1
                            # Now we want to save information about error's location
                            poses,muts=getErrors(primersR1_5[item[0]][1:-1],item[1][3])
                            for p in poses:
                                if p not in primersErrorsPos[0].keys():
                                    primersErrorsPos[0][p]=1
                                else:
                                    primersErrorsPos[0][p]+=1
                            for m in muts:
                                if m not in primersErrorsType[0].keys():
                                    primersErrorsType[0][m]=1
                                else:
                                    primersErrorsType[0][m]+=1
                        # Else it's a sequencing error
                        else:
                            primersErrors[item[0]][1][2]+=1
                    except IndexError:
                        print('IndexError!',a)
                        print(item)
                        exit(0)
                # Else we just save it as sequencing error
                else:
                    primersErrors[item[0]][0][2]+=1
            primersStatistics.write('Primer\tTotal_number_of_reads\tNumber_without_any_errors\t'
                                    'Number_with_sequencing_errors\tNumber_with_synthesis_errors\n')
            for key,item in primersErrors.items():
                item[0]=list(map(str,item[0]))
                item[1]=list(map(str,item[1]))
                primersStatistics.write(str(key+1)+'F\t'+'\t'.join(item[0])+'\n')
                primersStatistics.write(str(key+1)+'R\t'+'\t'.join(item[1])+'\n')
            primersStatistics.close()

            primersStatisticsPos.write('\t'.join(['Position_in_primer','Number_of_mutations_in_R1','Number_of_mutations_in_R2'])+'\n')
            for pos in range(1,max(max(primersErrorsPos[0].keys()),max(primersErrorsPos[1].keys()))+1):
                line=[str(pos)]
                if pos in primersErrorsPos[0].keys():
                    line.append(str(primersErrorsPos[0][pos]))
                else:
                    line.append('N/A')
                if pos in primersErrorsPos[1].keys():
                    line.append(str(primersErrorsPos[1][pos]))
                else:
                    line.append('N/A')
                primersStatisticsPos.write('\t'.join(line)+'\n')
            primersStatisticsPos.close()

            primersStatisticsType.write('\t'.join(['Error_type','Number_of_mutations_in_R1','Number_of_mutations_in_R2'])+'\n')
            for key in set(list(primersErrorsType[0].keys())+list(primersErrorsType[1].keys())):
                line=[key]
                if key in primersErrorsType[0].keys():
                    line.append(str(primersErrorsType[0][key]))
                else:
                    line.append('N/A')
                if key in primersErrorsType[1].keys():
                    line.append(str(primersErrorsType[1][key]))
                else:
                    line.append('N/A')
                primersStatisticsType.write('\t'.join(line)+'\n')
            primersStatisticsType.close()
        if idimer:
            idimerFile.write('Primer-dimer\tNumber of read pairs\n')
            for key,item in sorted(primerDimers.items(),key=itemgetter(1),reverse=True):
                idimerFile.write(key+'\t'+str(item)+'\n')
            idimerFile.close()

        trimmedReadsR1.close()
        untrimmedReadsR1.close()
        if args.trimmedReadsR2:
            trimmedReadsR2.close()
            untrimmedReadsR2.close()
    elif bamFile:
        print('Reading file with coordinates of amplicons...')
        coordToPrimerNum={}
        amplCoords=[]
        hardClipping=args.hardClipping
        file=open(coordsFile)
        primerNum=0
        for string in file:
            cols=string.replace('\n','').split('\t')
            if 'chr' not in cols[0]:
                cols[0]='chr'+cols[0]
            if cols[0] not in coordToPrimerNum.keys():
                coordToPrimerNum[cols[0]]={}
            for i in range(int(cols[1]),int(cols[2])+1):
                if i not in coordToPrimerNum[cols[0]].keys():
                    coordToPrimerNum[cols[0]][i]=[primerNum]
                else:
                    coordToPrimerNum[cols[0]][i].append(primerNum)
            amplCoords.append([cols[0],int(cols[1]),int(cols[2])])
            primerNum+=1
        file.close()
        print('Reading input BAM-file...')
        if bamFile[-3:]!='.gz':
            bam=pysam.AlignmentFile(bamFile)
            outBam=pysam.AlignmentFile(args.outBamFile,'wb',template=bam)
            if args.outUntrimmedBamFile:
                outBam2=pysam.AlignmentFile(args.outUntrimmedBamFile,'wb',template=bam)
            reads=bam.fetch()
        else:
            print('ERROR! cutPrimers does not support gzipped BAM-files:',bamFile)
            exit(0)
        p=ThreadPool(threads)
        # Cutting primers and writing result immediately
        print('Trimming primers from reads of BAM-file...')
        results=[]
        # Dictionary that stores positions of each read by read name and its flag value
        matePositions={}
        for r in reads:
            results.append(p.apply_async(trimPrimersInBam,args=(r,coordToPrimerNum,amplCoords,maxPrimerLen,primerLocBuf,errNumber,primersR1_5,primersR1_3,primersR2_5,primersR2_3,
                                                                primerR1_5_hashes,primerR1_5_hashLens,primerR2_5_hashes,primerR2_5_hashLens,
                                                                primersFileR1_3,primersFileR2_5,primersFileR2_3,primer3absent,minPrimer3Len,minReadLen,hardClipping)))
        doneWork=0
        allWork=len(results)
        if allWork!=0:
            showPercWork(0,allWork)
        newReads=[]
        for result in results:
            res=result.get()
            if res[0]:
                newReads.append(res[1])
                if res[1].is_read1:
                    matePositions[res[1].query_name+'_1']=[res[1].reference_id,res[1].pos]
                elif res[1].is_read2:
                    matePositions[res[1].query_name+'_2']=[res[1].reference_id,res[1].pos]
                else:
                    matePositions[res[1].query_name+'_0']=[res[1].reference_id,res[1].pos]
            elif args.outUntrimmedBamFile:
                try:
                    outBam2.write(res[1])
                except TypeError:
                    print('ERRROR!',res)
                    exit(0)
                if res[1].is_read1:
                    matePositions[res[1].query_name+'_1']=[None,None]
                elif res[1].is_read2:
                    matePositions[res[1].query_name+'_2']=[None,None]
                else:
                    matePositions[res[1].query_name+'_0']=[None,None]
            doneWork+=1
            showPercWork(doneWork,allWork)
        print()
        p.close()
        p.join()
        bam.close()
        if len(newReads)==0:
            print('WARNING! No reads left after cutting primer sequences! Possibly, your reads do not have primer sequences on the 3\'-ends. If so, use parameter, -primer3')
            os.rename(args.outBamFile,args.outBamFile[:-4]+'.sorted.bam')
        else:
            print('Changing coordinates of mates...')
            p=ThreadPool(threads)
            results=[]
            for r in newReads:
                results.append(p.apply_async(changeMateCoordinates,args=(r,matePositions)))
            doneWork=0
            allWork=len(results)
            showPercWork(0,allWork)
            for result in results:
                res=result.get()
                outBam.write(res)
                doneWork+=1
                showPercWork(doneWork,allWork)
            print()
            p.close()
            p.join()
            outBam.close()
            print('Sorting output BAM-file...')
            pysam.sort('-o',args.outBamFile[:-4]+'.sorted.bam',args.outBamFile)
            if args.outUntrimmedBamFile:
                outBam2.close()
                pysam.sort('-o',args.outUntrimmedBamFile[:-4]+'.sorted.bam',args.outUntrimmedBamFile)
            print('Indexing output BAM-file...')
            pysam.index(args.outBamFile[:-4]+'.sorted.bam')
            if args.outUntrimmedBamFile:
                pysam.index(args.outUntrimmedBamFile[:-4]+'.sorted.bam')
        




        

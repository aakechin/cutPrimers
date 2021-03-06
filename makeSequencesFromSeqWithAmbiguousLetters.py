# This script makes all possible pairs of primers from sequences with abmiguous letters

import argparse
from Bio import Seq

if __name__ == "__main__":
    # Section of reading arguments
    par=argparse.ArgumentParser(description='This script makes all possible sequences from sequence with abmiguous letters')
    par.add_argument('--inputFastaFile','-in',dest='inFastaFile',type=str,help='fasta-file with sequences of primers for all amplicons. First, you should write forward primer, then reverce one for one amplicon. Second, in the same order for the 2nd amplicon. And then in the same order for next amplicons',required=True)
    par.add_argument('--output','-out',dest='outFile',type=str,help='basis for an output-file',required=False)
    args=par.parse_args()

ambNucs={'W':['A','T'],'K':['G','T'],'M':['A','C'],'S':['C','G'],'R':['A','G'],'Y':['C','T'],
        'B':['C','G','T'],'V':['A','C','G'],'D':['A','G','T'],'H':['A','C','T']}

file=open(args.inFastaFile)
seqNames=[]
seqs=[]
for string in file:
    if '>' in string:
        seqNames.append(string.replace('>','').replace('\n',''))
    else:
        seqs.append(string.replace('\n',''))
file.close()

newSeqsF=[]
newSeqsR=[]
for seqNum,seq in enumerate(seqs):
    ambNucPos=[]
    varNum=1
    if seqNum%2==0:
        newSeqsF.append([])
        for i,c in enumerate(seq):
            if c in ambNucs.keys():
                ambNucPos.append(i)
                varNum*=len(ambNucs[c])
        for i in range(varNum):
            stats=str(bin(i))[2:].zfill(len(ambNucPos))
            newSeq=''
            for j,c in enumerate(seq):
                if c in ambNucs.keys():
                    newSeq+=ambNucs[c][int(stats[ambNucPos.index(j)])]
                else:
                    newSeq+=c
            newSeqsF[-1].append(newSeq)
    else:
        newSeqsR.append([])
        for i,c in enumerate(seq):
            if c in ambNucs.keys():
                ambNucPos.append(i)
                varNum*=len(ambNucs[c])
        for i in range(varNum):
            stats=str(bin(i))[2:].zfill(len(ambNucPos))
            newSeq=''
            for j,c in enumerate(seq):
                if c in ambNucs.keys():
                    newSeq+=ambNucs[c][int(stats[ambNucPos.index(j)])]
                else:
                    newSeq+=c
            newSeqsR[-1].append(newSeq)

rFileR1_5=open(args.outFile+'_R1_5.fa','w')
rFileR1_3=open(args.outFile+'_R1_3.fa','w')
rFileR2_5=open(args.outFile+'_R2_5.fa','w')
rFileR2_3=open(args.outFile+'_R2_3.fa','w')
for k,(seqsF,seqsR) in enumerate(zip(newSeqsF,newSeqsR)):
    for i,seqF in enumerate(seqsF):
        for j,seqR in enumerate(seqsR):
            rFileR1_5.write('>'+seqNames[2*k+1]+'_'+str(j+1)+'\n'+seqR+'\n')
            rFileR2_5.write('>'+seqNames[2*k]+'_'+str(i+1)+'\n'+seqF+'\n')
            rFileR2_3.write('>'+seqNames[2*k+1]+'_'+str(j+1)+'_reverse_complement\n'+str(Seq.Seq(seqR).reverse_complement())+'\n')
            rFileR1_3.write('>'+seqNames[2*k]+'_'+str(i+1)+'_reverse_complement\n'+str(Seq.Seq(seqF).reverse_complement())+'\n')
rFileR1_5.close()
rFileR1_3.close()
rFileR2_5.close()
rFileR2_3.close()

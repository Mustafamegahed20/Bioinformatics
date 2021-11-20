import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
import numpy as np
import re

def Excludeintrons(inputFileName, intronsFileName):
    dnaFile = open(inputFileName)
    dna = dnaFile.read()
    dnaFile.close()
    print('Before transcription:', dna)


    intronsPosFile = open(intronsFileName)
    intronsPos = intronsPosFile.read()
    intronsPosFile.close()

   
    pos = intronsPos.split("\n")

    newdna = dna

    for p in pos:
        if p != "":
            startEndPos = p.split(",")
            print('positions:', p)
            newdna = newdna.replace((dna[int(startEndPos[0]):int(startEndPos[1])+ 1]), '')

    return newdna


returned = Excludeintrons("genomic_dna.txt", "introns.txt")
print('After transcription:', returned)



def getFragmentSize(seq):

    s= re.search(r"GAATTC",seq)
    fregment1 = 1 + s.start()

    fregment2 = len(seq) - fregment1

    return fregment1, fregment2


Fregment1, Fregment2 = getFragmentSize("ACTGATCGATTACGTATAGTAGAATTCTATCATACATATATATCGATGCGTTCAT")
print(Fregment1)
print(Fregment2)

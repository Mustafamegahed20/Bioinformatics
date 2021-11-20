import re
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
import numpy as np


def getFragmentSize(seq):

    fregment1 = 1 + seq.find("GAATTC")

    fregment2 = len(seq) - fregment1

    return fregment1, fregment2


Fregment1, Fregment2 = getFragmentSize("ACTGATCGATTACGTATAGTAGAATTCTATCATACATATATATCGATGCGTTCAT")
print(Fregment1)
print(Fregment2)
